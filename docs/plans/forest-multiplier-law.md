# forest-multiplier-law: what a forest's coefficient is, held or drawn, and what a numeric multiplier's sd is the size of

Status: PLANNED (dec-B253 and dec-B246 as revised; dec-B260 to dec-B265; dec-B269, dec-B271, dec-B272, dec-B275,
dec-B276; dec-A171). Follows [forest-kind-by-class.md](forest-kind-by-class.md),
[forest-sd-unit.md](forest-sd-unit.md), [forest-defaults-by-kind.md](forest-defaults-by-kind.md) and push 3 of
[written-surface.md](written-surface.md), none of which has landed.

One plan, four pushes. The engine's change of law and the R surface that states it are not planned apart: at
every tip the help must say what the engine does and the reader must report the number the engine uses, and
a plan for the engine alone would leave a tip that does neither.

agent: push 1 (statements): opus implementer for the engine and the bridge, sonnet for the respelled pins,
opus reviewer. Pushes 2 and 3 (the law): opus implementer for the engine, the bridge and the R code, sonnet
for the respelled tests, the benchmark scripts and the help once the code is fixed, opus reviewer told to
refute. Push 4 (the printed block): sonnet implementer, opus reviewer. The reason for opus on 1 to 3: every
slip is a prior off by a factor and no message. The scale of a column can be taken from the wrong rows or
taken twice on any of nine paths, a held value can go by position on one path and by kind on another, and a
number per column can become one number anywhere between R and the draw; a stand-in for the slice, built
with care, had three such faults (Context).
rng: stated per push and per call sequence, as [RNG classes and their gates](README.md#rng-classes-and-their-gates)
defines the classes.
- Push 1: NEUTRAL for every call both builds accept. A sampler saved before it is refused by name when it is
  next used.
- Push 2: POSTERIOR-CHANGING for a model that holds a coefficient, in four sequences: a held two-level factor
  (the forest's own sd is s where it was s / 0.674); a held factor that is not the second forest (held at 0, 1
  where it was held at 1, 1); a held forest with no basis that states an sd, or is left at its default under
  a gaussian response (its sd is s where it was the response's scale whatever s; under probit and logistic
  the default is the same number and the fit keeps its bits); a held forest with no basis that is not the
  first forest (held at 1 where it was held at 0). NEUTRAL for every model whose coefficients are all drawn,
  numeric multipliers included. Refused where accepted: `amplitude.prior.variance`, a held factor of three or
  more levels, a held numeric basis of two or more columns.
- Push 3: POSTERIOR-CHANGING for every model with a numeric multiplier, drawn, at the default sd or a stated
  one. NEUTRAL for every model whose forests have no basis or a factor, held or drawn, stated or not.
  Accepted where refused: a held coefficient on one numeric column. Refused where accepted: a numeric basis
  with a constant column and no sd; a swap that changes the width of a numeric forest at its default sd.
- Push 4: no draw moves.
- One recorded scenario moves, in push 2: `glue_toggle` of the BCF equivalence baseline, which holds the
  treatment forest's coefficients. No other scenario of the three baselines holds a coefficient or multiplies
  a number, and the four snapshot files fit no model of several forests (searched; the BCF compare and the
  gaussian baseline's one scenario of several forests were run on a stand-in: Context).
window: pre-release, after forest-kind-by-class and before the control-migration arc moves the forests'
record to the model. Serial with any other work in [`ForestSpec`](../../src/bartcore/combiner.hpp),
[`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp), the K-forest constructor of
[`Chain`](../../src/bartcore/chain.hpp), [`Chain::setForestBasis`](../../src/bartcore/chain.hpp),
[`Chain::setForestMapSd`](../../src/bartcore/chain.hpp),
[`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp), [`parseData`](../../src/R_interface_bartcore.cpp),
[`bartcore_getLeafPrior`](../../src/R_interface_bartcore.cpp), [`forestParams`](../../R/model.R), the
multi-forest block of [`resolveSamplerSpec`](../../R/spec.R), [`reportLeafPrior`](../../R/dbarts.R) or the
sampler's `initialize`; so serial with [leaf-conversions.md](leaf-conversions.md) and
[cross-family-state-install.md](cross-family-state-install.md), in either order (What waits on what).
bartCause's edit lands the day push 2 does.
budget: ~7800 lines changed, upper figure 12500. By push: 1 ~1250 (upper 2050), 2 ~1950 (upper 3150),
3 ~4000 (upper 6350), 4 ~600 (upper 950). By layer, with the upper figure: engine 770 (1260); bridge 350
(580); tests/cpp 1140 (1850); R 1030 (1660); tinytest 2430 (3900); benchmarks 900 (1440); help 520 (840);
design note, architecture, TODO and the two indexes 630 (980). The design estimated 3600 and planned for
6200. This figure is counted, not scaled: the stand-in's engine and bridge are 650 changed lines with no
comment in the house style, its run of the suite shows 121 assertions of this slice's to repair and 29
creations to respell, and the exact gate's draft is 150 dense lines; written-surface's two landed pushes
ran at 1.8 and 2.7 times their plans.

## Goal

Each kind of forest has one law, and it is the law the help states. A forest's `sd` is the size of what it
contributes for one step of its multiplier: exactly where its coefficient is held, as the prior median where
it is drawn. A numeric column is used as it was given, never centred and never divided; its coefficient is
normal; a stated sd is per unit of the column; and with no sd stated the default is sqrt(2 / K) of the
response's scale for one standard deviation of each column, that standard deviation taken once, when the
sampler is created, and kept with the data. A held coefficient has the value its forest's kind gives it,
wherever the forest stands, and the three shapes that have such a value are the three accepted.
`amplitude.prior.variance` is not an argument. The reader, `extract` and `print` give one number per column
of a basis, and each is the number the engine draws under.

## Context

Measured at the tip (d02fe72c; written-surface pushes 1 and 2 landed; push 3, forest-defaults-by-kind,
forest-sd-unit and forest-kind-by-class not) on the shipped build, R 4.6.1. The spellings are the tip's, a
tilde on a basis written as code. L is the response's scale as the engine holds it: sd(y - offset) under
gaussian, 1 under probit, pi / sqrt(3) under logistic; forest-sd-unit makes it the unit an sd is divided by
and changes nothing else below. K is the number of forests. "Size" is the prior median of the absolute
multiplier times the prior standard deviation of the forest's own fit, from the engine's prior draws under
a flat likelihood (4000 sweeps a row, Monte Carlo error about 0.02); a held multiplier is exact. Fixture:
300 rows; x1, x2; a 0/1 z; factors of 3 and 8 levels; w with standard deviation 3.15; two columns with
standard deviations 2.0 and 0.48; an age with mean 50.

- The suite. 231 files, of which 4 exit off a reference build and 3 ask for more than two threads; the other
  224 give 16065 results, none failing, 96 seconds in one process.
- The law today. s is the stated or default number, c the median nonzero row norm of the basis.

  | forest | coefficient | the forest's own sd | what was measured, over s L |
  |---|---|---|---|
  | no basis, drawn | N(0, v), v inverse-gamma: its size half-Cauchy, median s (measured 2.05 at s = 2; q90 / q50 6.6, the law's 6.3) | L (0.999) | a F: 1.02 |
  | no basis, held | 1 | L, whatever s (0.995) | F: 1.42 at s = 0.7; 0.50 at the gaussian default 2; 0.99 at the latent default 1 |
  | factor of 2, 3 or 8 levels, drawn | N(0, 1/2) each (sd 0.71) | s L / 0.674 (1.48 at s = 1) | (a_k - a_l) F: 0.99 to 1.01 |
  | two-level factor, held, second forest | (0, 1) | s L / 0.674 | (a_2 - a_1) F: 1.49 |
  | one numeric column, drawn | N(0, 1/2) | s L / (0.674 c) | a F per unit: 0.72 / c; per sd of the column 1.12 for w and 0.24 for the age |
  | two numeric columns, drawn | N(0, 1/2) each | s L / (0.674 c), one c for both | each a_j F per unit: 0.69 / c |
  | two numeric columns, held | (0, 1) | s L / (0.674 c) | the second column's term: 1.42 / c; the first column's is zero |
  | one numeric column, held | | | refused (dec-A171) |

  The defaults: 2 for a forest with no basis under gaussian and 1 under probit and logistic; sqrt(2 / K) for
  one with a basis, a factor and numbers alike. The same figures hold under probit and logistic (nine rows
  measured, each the gaussian row's to Monte Carlo error).
- `amplitude.prior.variance = v`, on a forest with a basis: the coefficients are N(0, v) and the forest's own
  sd does not move, so the size is sqrt(2 v) times the table's: 1.95 at v = 2 and 0.52 at v = 1/8 on a factor
  (2.04 with `sd = 0.7` beside it), 1.03 / c at v = 1 on a number. Nothing on a held forest. Refused on a
  forest with no basis. Stated on 31 lines of 6 test files, in 5 benchmark scripts and by bartCause.
- A held value goes by position. [`AmplitudeState`](../../src/bartcore/combiner.hpp) starts at 1 for the
  first forest and (0, 1) for the second, and [`rebuildAmplitudeLayout`](../../src/bartcore/combiner.hpp)
  carries those by position and fills every other coefficient with 1. Read off the engine: a factor of three
  levels held as the second forest is (0, 1, 1), its third level's term the second's; of eight levels (0, 1,
  1, 1, 1, 1, 1, 1); a two-level factor held as the third forest, or as the first where no forest is without
  a basis, is (1, 1): no contrast, a function added on every row (size of the difference 0.000). A drawn
  coefficient starts at the same values. In a tests/cpp fixture a forest with no basis held as the SECOND
  forest is held at 0, with one column or two before it: R cannot write that model at the tip (a forest past
  the first needs a basis) and can once forest-defaults-by-kind lands.
- The row norm. c is taken over the rows the sampler holds, after `subset` (1.976 for rows 1 to 200, R's
  median of those rows), over the rows where the basis is not zero, unweighted; rows of weight zero and rows
  masked later count. It is taken again at every swap: `$setForestBasis(2, 3 * w)` divides the forest's own sd
  by 3 (2.493 to 0.831) and leaves the coefficient; a swap to two columns enters the new coefficient at 1.
  `copy()`, a reload and `new("dbartsSampler", control, model, data)` take it from the column then in force,
  so after a swap to w / 10 all three have ten times the creation's sd. Nothing of it is recorded.
- Other paths. `setData` is refused on every sampler whose forests carry coefficients, and
  `setResponse(updateScale = TRUE)` too; `setResponse(updateScale = FALSE)`, `setOffset`, `setWeights` and
  `setPredictor` move no forest's sd. A control taken from a sampler and given to a fit of a response twenty
  times as wide carries no anchor: the new sampler's is its own (67.71 against 3.386). A warm start refuses
  several forests, and a model of several forests refuses test predictors, so no test row can reach a basis.
- The prior-only path. `sampleTreesFromPrior` with `sampleLeafParametersFromPrior` draws each forest's own
  fit at the sd the reader reports (three forests over 3000 draws: 0.997, 1.027, 0.987 of `k.scale`) and
  leaves every coefficient where it was. [`samplePriorPredictive`](../../R/dbarts.R) on a sampler of several
  forests stops inside `predict` ("no off-sample basis"), with a formula term too.
- A state. Its glue block is K, the widths, every coefficient, held ones included, and one variance per
  forest ([`readAmplitudeGlue`](../../src/R_interface_bartcore.cpp)). An install takes a forest's
  coefficients where the recipient draws them and its variance where the recipient's size is half-Cauchy
  ([`restoreGlue`](../../src/bartcore/combiner.hpp)); a held block in the state is passed by, the rule
  ["fixed amplitudes"](../../inst/tinytest/test-state-not-model.R) pins. So a state from a drawn sampler goes
  into one that holds (0, 1), the trees arriving and the coefficients staying (0, 1), "exact" `TRUE`. A state
  holds no sd, no scale, no row norm and no kind.
- The engine. [`ForestSpec`](../../src/bartcore/combiner.hpp) carries four derived numbers per forest
  (`nodeScaleFactor`, `nodeScaleDivisor`, `amplitudePriorVariance`, `amplitudePriorScale`) and the hold flag;
  R derives them ([`forestParams`](../../R/model.R)) and the bridge passes eight numbers a forest
  ([`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp)). The K-forest constructor of
  [`Chain`](../../src/bartcore/chain.hpp) divides by [`basisRowNorm`](../../src/bartcore/chain.hpp)
  ([`mapLeafScale`](../../src/bartcore/chain.hpp)). [`ForestAmplitudePrior`](../../src/bartcore/combiner.hpp)
  holds one variance for a forest's whole block, and
  [`drawForestAmplitude`](../../src/bartcore/combiner.hpp) seeds every coordinate of the block's precision
  with it: the engine cannot give two columns of one forest two scales.
- The shipped header. [`DBARTS_C_API_LIST`](../../inst/include/dbarts/dbarts.h) creates no sampler from a
  forest's statement and reads no forest's spread, scale or hold; a per-draw callback sees the coefficients
  as `glue`. Nothing in this slice changes a signature there.
- The reader. `$getLeafPrior(f)` on a forest of several gives `leaf.prior` (`forest(sd = s)`),
  `prior.sd.of` ("amplitude scale" or "forest total"), `k.scale`, and `amplitude.prior.scale` or
  `amplitude.prior.variance`, `leaf.scale.factor`, `leaf.scale.divisor`, `basis.row.norm`; with no forest, an
  unnamed list. `extract(type = "leaf.prior.sd")` is a vector `forest1`, `forest2`, each the forest's own sd.
  `print` of a fit, `show` of a sampler and the verbose summary say nothing of any forest's sd.
- A stand-in for the slice: the forest-prior-args design's prototype of the law, ported to the tip (the
  engine, the bridge and R in the tip's spelling, an sd still in units of L; it also carries the
  kind-by-class rule). On it:
  - The law. 28 rows of the table under "The rule", every kind, held and drawn, default and stated, at the
    second, third and first position, under the three families: measured size over the rule's 0.967 to
    1.032. tests/cpp builds against it and passes, 350 lines ok.
  - What keeps its bits. 20 of 20 seeded fits of forests with no basis or a factor, coefficients drawn, are
    `identical()` on the tip and the stand-in (two and three forests, a factor third, 8 levels, stated sds,
    probit, logistic, weights and an offset, `subset`, `bart` with two chains, a logical basis, a swap, a
    write, a copy, a reload, `setResponse`, forest weights and a mask, no forest without a basis); a held
    forest with no basis at its default under probit too. 11 fits with a number or a hold differ. The BCF
    equivalence compare: 14 of 15 scenarios identical, `glue_toggle` not (111 of 209 summaries beyond z = 3,
    the largest 13.95); the one scenario of several forests in the gaussian baseline is bitwise.
  - The suite. As the stand-in stands 7 files stop at a refusal and 25 assertions fail. With its refusals
    logged and passed by, 15972 results run, one file stops and 132 fail in 9 files (test-bcf-family.R 76,
    test-forest-arguments.R 25, test-bcf-creation.R 13, test-calibration-midchain.R 6,
    test-multiforest-leaf-prior-writer.R 6, four files 6 between them): 72 pin the stored numbers by
    position, 28 the reader's entries and labels, 13 the row norm, the coefficient's variance and a start by
    position, 8 are dec-A171's refusals, 11 are the kind-by-class slice's. No seeded draw of a drawn forest
    with no basis or a factor moved. The refusals passed by: 29 creations on a constant
    column with no sd (test-bcf-family.R 27), one held numeric basis of two columns, and 20 swaps the
    kind-by-class slice respells.
  - Its faults, each a requirement here. A declaration that gives a forest of a data object another basis
    (`dbartsSpec()` over `sampler$data`) keeps the old column's recorded scale: sd 0.317 where the new column,
    a hundred times as wide, has 0.00317. A forest stated and then swapped to two columns reports sds of
    1e-172, and its `copy()` stops on the stale record. On the column 1e8 + 1e-8 k the engine's one-pass
    standard deviation is 1.02e-6 where R's is 8.68e-7.
- The exact gate, drafted and run against the prototype's own build. One binary predictor, one tree a forest,
  sigma held; the leaves integrate out and the only quadrature is over the coefficients. Eleven arms (a 0/1
  column; a column of standard deviation 0.107; two columns of 0.29 and 7.0; a column with a fifth of its
  rows at weight zero; default and stated; drawn and held) at 8 seeds of 100000 kept sweeps: every gap under
  0.0043, 67 seconds on two cores for the eleven. Against a sampler given a wrong law through its record, and
  against two engine mutants: a default per unit 0.06 to 3.2; both columns on the first one's scale 1.36; the
  two scales exchanged beyond 1; the scale taken over the rows of positive weight 0.18 and 0.47; over the
  nonzero rows 0.047, and a 0/1 column refused; a precision divided by the scale and not its square 0.79.
  Against an oracle written to a wrong law: 0.674 left on a held coefficient 0.03 at the default, 0.10 and
  0.21 at sds of 0.6 and 0.2 of the unit; a coefficient variance of 1/2, 0.11; a held coefficient at 0, 0.27;
  a centred column 0.51; a stated sd read per standard deviation 0.14 to 0.65; a held plain forest that
  ignores an sd of 0.3 of the unit, 0.033. Not seen: n for n - 1 in the standard deviation (0.004 at 40
  rows).
- What the BCF gates hold. [bcf-exact.R](../../benchmarks/R/bcf-exact.R),
  [bcf-exact-weak.R](../../benchmarks/R/bcf-exact-weak.R),
  [bcf-exact-restricted.R](../../benchmarks/R/bcf-exact-restricted.R) and
  [bcf-latent-exact.R](../../benchmarks/R/bcf-latent-exact.R) each hold a coefficient in some arm, write
  ["scaleTau"](../../benchmarks/R/bcf-exact.R) with the 0.674 whether it is held or not, and pass
  `amplitude.prior.variance = 0.5`; [sbc.R](../../benchmarks/R/sbc.R) has a
  ["fixedGlue"](../../benchmarks/R/sbc.R) arm. Run on the stand-in, in `quick`, with their oracles as they
  are: the three gaussian gates pass under the new law (largest gaps 0.0034, 0.0224 against 0.03, and
  0.0004), so they cannot see it; the latent gate fails in six of the seven arms that hold the treatment
  forest's coefficients (worst |z| 5.34 to 10.70; the seventh, three forests under probit, 3.97 and inside
  its bound) and passes in the two that hold only the prognostic one at its default (1.72, 2.66), in 216
  seconds.
- Consumers. bartCause's `bcf()` (dbarts-1.0, 6c1bff9) passes `amplitude.prior.variance = b.prior.variance`
  on its treatment forest and holds either forest with `amplitude = fixed()`; it reads
  `getLeafPrior(1L)$response.scale` and `$response.shift` and the coefficients as `glue`, and nothing else of
  this surface. Two of its test files build the same sampler by hand with `amplitude.prior.variance = 0.5`.
  stan4bart (bartcore, a9d081b), treatSens (dbarts-1.0, aecec71) and bairrtt (main, 3f57f61) declare one
  forest with `dbartsForests$forest(n.trees = )` and read no coefficient, no forest's sd and no leaf prior of
  a model of several forests.

## The rule

s is a forest's `sd`, in the response's units; U is the unit forest-sd-unit records (sd(y - offset) over the
rows kept under gaussian, 1 under probit, pi / sqrt(3) under logistic); d is sqrt(2 / K). Term f at row i is
(B_f(i, .) . a_f) F_f(x_i): the basis row times the forest's coefficients, times the forest's own fit. The
kind of a forest is what forest-kind-by-class records: nothing, the levels of a factor, or numeric columns.

| the forest multiplies | coefficients | drawn | held, `amplitude = fixed()` |
|---|---|---|---|
| nothing | one | F has sd U; the size of a is half-Cauchy with median s / U: a F has median size s | a is 1; F has sd s, exactly |
| the q levels of a factor | one a level | each a_l is N(0, 1/2); F has sd s / 0.674: (a_k - a_l) F has median size s | q = 2 only: a is (0, 1); F has sd s, exactly, the second level against the first |
| numeric columns w_1 to w_q | one a column | a_j is N(0, (s_j / s_1)^2), each alone; F has sd s_1 / 0.674: a_j F has median size s_j, per unit of w_j | q = 1 only: a is 1; F has sd s_1, exactly, per unit of w |

1. Stated. `sd = s` is one number. On a numeric basis every s_j is s, per unit of each column as it was
   given, and nothing is read from the data.
2. Not stated. No basis: 2 U under gaussian, U under probit and logistic. A factor: d U. Numbers:
   s_j = d U / sd(w_j), each column by its own standard deviation. Held or drawn alike.
3. The scale, sd(w_j): the sample standard deviation of the column, n - 1 in the divisor, unweighted, over
   every row the sampler holds when it is first created. Those are the rows left by `subset` and the
   na.action; rows of weight zero count, and no later mask enters. The engine takes it, once, and
   R records it with the data; every later construction is handed the record and takes nothing again. A
   constant column has no scale: where a default would need one it is refused by name.
4. The column is never centred and never divided: the engine multiplies by it as given, at the fit, after a
   swap and at new rows. A default fit on w and on c w is the same fit.
5. A held value goes by the forest's kind and never by where the forest stands. The three held shapes of the
   table are the three accepted; every other is refused by name, in R and again in the engine.
6. The variance of a coefficient is not an argument: `amplitude.prior.variance` goes.
7. A swap moves the basis and nothing else: no scale, no sd in force, no coefficient, no tree. A numeric
   forest at its default sd keeps its width: a column added would have no default.
8. A write, `$setLeafPrior(forests = )`, states: s for every column, and the forest is a stated one from
   then on, also where the number written is the one in force, in which case no bit of the prior moves.
9. Where a drawn coefficient starts is not part of the law. A forest with no basis or a factor starts where
   it starts today, which is what keeps its draws; a numeric coefficient starts at s_j / s_1.

Before and after, per path.

| path | before | after |
|---|---|---|
| creation at every door: a `forest()` term at `bart()` and `dbarts()`; a `forests` list on a formula, on a matrix, in `dbartsSpec()`; `dbartsData(bases = )` | the table of Context | the rule; a numeric forest that states no sd has its scale taken and recorded |
| creation from a data object that carries a scale (`sampler$data`) | | used, for a forest that states no sd; its width and sign checked there and nowhere else |
| a declaration that gives a forest of such an object another basis | | that forest's record is dropped and taken from the new column |
| `$setForestBasis` | the row norm is taken again: the forest's own sd follows the new column | nothing but the basis moves; a default numeric forest refuses another width |
| `$setLeafPrior(forests = )` | the number goes to one of two channels | the forest is stated at s for every column |
| `copy()`, a reload, `new("dbartsSampler", control, model, data)` | built from the column then in force | built from the record: the creation's prior, bit for bit |
| `setState`, a state in a reload | carries no prior | unchanged; a held block in the state is passed by, as today |
| `predict` at new rows, `fitted`, `extract` of draws | coefficients times the basis as given | unchanged |
| `sampleTreesFromPrior`, `sampleLeafParametersFromPrior` | each forest's own sd; coefficients left | unchanged, at the rule's own sd |
| `samplePriorPredictive` | stops inside `predict` | refused by name |
| `setData`; a warm start | refused | refused |
| `setResponse`, `setOffset`, `setWeights`, `setPredictor`, `setActiveRows` | move no forest's sd | unchanged |
| a control or a `forests` list carried to another fit | the forests' record is cleared | unchanged; a scale is on the data and never on the control |
| a sampler or a fit saved before push 1 | | refused by name when it is next used |
| a state stored before | | installs as any state; it holds nothing of the law |

## The reader, `extract` and the printed block

One forest's entry from `$getLeafPrior(f)` after push 3, on a model of several forests. The figures are for a
response with sd 1.452, a dose with sd 0.0204 and an age with sd 17.55; in the table each forest stands
beside one forest with no basis (K = 2), and the list and the block below are one model of four (K = 4).
These are the returns the tests pin, each literal computed in the test from its fixture.

| forest | `leaf.prior` | `sd` | `sd.stated` | `multiplier` | `amplitude` | `prior.sd.of` |
|---|---|---|---|---|---|---|
| no basis, default | `forest()` | `2.904` | `FALSE` | "none" | absent | "forest (prior median)" |
| one numeric column, default | `forest()` | `c(dose = 71.18)` | `FALSE` | "numeric" | absent | "forest, per unit of basis (prior median)" |
| the same, `sd = 30` | `forest(sd = 30)` | `c(dose = 30)` | `TRUE` | "numeric" | absent | the same |
| two numeric columns, default | `forest()` | `c(dose = 71.18, age = 0.08274)` | `FALSE` | "numeric" | absent | the same |
| the same, `sd = 30` | `forest(sd = 30)` | `c(dose = 30, age = 30)` | `TRUE` | "numeric" | absent | the same |
| factor, default | `forest()` | `1.452` | `FALSE` | "factor" | absent | "level difference (prior median)" |
| two-level factor, held | `forest()` | `1.452` | `FALSE` | "factor" | `fixed()` | "level difference" |
| one numeric column, held, `sd = 30` | `forest(sd = 30)` | `c(dose = 30)` | `TRUE` | "numeric" | `fixed()` | "forest, per unit of basis" |

- `leaf.prior` is the statement: what was stated, and the default as the default (dec-B269). It goes back
  into `$setLeafPrior(forests = )` as it is, and a default written back states nothing. `sd` is the number
  in force, one for each column of a numeric basis, named as the column is where push 3 of written-surface
  names it and unnamed in column order otherwise; one unnamed number for a forest with no basis or a factor.
- `basis.scale`: on a numeric forest whose data carries one, the recorded standard deviations, named as `sd`.
- Kept as they are: `leaf.model`, `prior.mean`, `k.scale` (the forest's own sd), `response.scale`,
  `response.shift`. Gone: `amplitude.prior.variance` (push 2); `amplitude.prior.scale`, `leaf.scale.factor`,
  `leaf.scale.divisor`, `basis.row.norm` (push 3).
- With no forest the reader's list is named by label
  ([The label rule](written-surface.md#the-label-rule)).
- `extract(fit, type = "leaf.prior.sd")` on a fit of several forests is the list of those `sd` entries by
  forest, named by label (dec-B275):
  `list(forest1 = 2.904, dose = c(dose = 50.33), "dose + age" = c(dose = 50.33, age = 0.0585), "factor(z)" = 1.027)`.
  `forest =` selects entries of the list and the result is a list whatever is selected. On one forest it is
  what it is today.
- The printed block, in `print` of a fit, `show` of a sampler and the verbose summary (push 4):

      forests, sd in the response's units (its sd 1.452):
        1 forest1: no basis; sd 2.904 (default)
        2 dose: numeric basis, mean 0.01994, sd 0.0204; sd 50.33 per unit, 1.027 per sd (default)
        3 dose + age: numeric basis, 2 columns (default)
            dose  mean 0.01994, sd 0.0204; sd 50.33 per unit, 1.027 per sd
            age   mean 50.88, sd 17.55; sd 0.0585 per unit, 1.027 per sd
        4 factor(z): factor basis, 2 levels; sd 1.027 between two levels (default); coefficient held

  A stated sd reads `sd 30 per unit, 0.612 per sd (stated)`. Under probit the header is `forests, sd on the
  latent scale (error sd 1):` and under logistic `(error sd 1.814)`. A column's mean and sd are those of the
  column in force; where a default was taken from another column the line ends `(default, taken from a
  column of sd 0.0204)`. A mean is printed through `zapsmall` against the sd. This is the one sign a 0/1
  number (dec-B263: at 20 percent treated beside one other forest, `mean 0.1967, sd 0.3981; sd 3.647 per
  unit` where the factor's line says 1.452), a column far from zero (dec-B264: `mean 50.88, sd 17.55`) and a
  stated sd in the wrong units get: nothing warns.

## Refused forms, with their texts

Base R's style, as in written-surface. `<f>` is `forest 2`, or `forest 2 ("dose")` where the forest has a
label. The push that adds each is in brackets.

    [1] this sampler was saved by a build of dbarts that recorded its forests another way (8 numbers a forest, where 4 are read), and it cannot be run under the prior it had; create it again
    [2] <f> cannot hold its coefficients (amplitude = fixed()): its basis has 3 levels, and no held value sets every difference between levels; hold a two-level factor, or let them be drawn
    [2] <f> cannot hold its coefficients (amplitude = fixed()) on a numeric basis of 2 columns: no held value is defined for several columns; let them be drawn
    [2] <f>: amplitude = fixed() on a basis of one numeric column is not supported yet. Let the coefficient be drawn; a column of two values can be held if it is written as a factor
    [2] 'amplitude = fixed(2)': a held coefficient is 1 for a forest with no basis, and 0 and 1 for the two levels of a factor; fixed() takes no other value. Write fixed(), and state the forest's size with 'sd'
    [3] <f> cannot hold its coefficients (amplitude = fixed()) on a numeric basis of 2 columns: a held forest multiplies one column; give it one, or let them be drawn
    [3] 'amplitude = fixed(2)': a held coefficient is 1 for a forest with no basis and for one numeric column, and 0 and 1 for the two levels of a factor; fixed() takes no other value. Write fixed(), and state the forest's size with 'sd'
    [3] <f> states no sd and column 1 of its basis is constant, so there is no standard deviation to take a default from; state one with sd =
    [3] <f> has a default sd for each of its 1 basis columns; a basis of 2 columns needs a stated sd: state one with $setLeafPrior
    [3] <f> states no sd and the standard deviations of its basis columns are too far apart to give each a default (1e-100 and 1e+100); rescale a column in 'basis', or state an sd
    [3] samplePriorPredictive does not support a sampler of several forests: their coefficients are not drawn from the prior here
    [3] 'basis.scale' entry 2 must hold one number for each column of that forest's basis: it has 1 and the basis 2 columns            ("must be finite and positive")
    [3] cannot extract 'leaf.prior.sd': this fit of several forests was saved before fits recorded an sd for each column of a basis; fit it again

An argument that does not exist gets R's own error: `unused argument (amplitude.prior.variance = 0.5)` from
push 2, and `updateBasisScale =` on `$setForestBasis` as today. Behind these, from the engine and the bridge:
"a held coefficient block must be a forest's with no basis, a two-level block or one numeric column"; "a
forest with no basis takes none"; "a basis column is constant, so no default sd can be taken from its
standard deviation"; "an sd for each column needs a numeric basis of that many columns"; "the sampler refuses
this basis for the forest", where today the bridge drops the engine's answer. Gone: dec-A171's refusal of one
numeric column, with push 3 (its first form goes with push 2);
"'amplitude.prior.variance' is the prior on a basis forest's amplitudes", with push 2. Kept as they are: every
text of the kind-by-class slice, whose text for a held forest's width stands.

## The help's `sd` item, as push 3 leaves it

    sd: The size of what this forest contributes, in the response's units: the units of y for a
    continuous response, and of the latent index under "probit" and "logistic", where the link's own
    error has standard deviation 1 and pi / sqrt(3). One number. What it is the size of depends on what
    multiplies the forest.

    No basis. The forest enters as a f(x), and sd is the size of a f(x).
    A factor (a factor, a character vector or a logical vector). The forest enters with one coefficient
    for each level, and sd is the size of the difference between two levels: the treatment effect where
    the levels are control and treated.
    Numbers w. A column that is an ordinary covariate should be centred: basis = scale(w), or
    basis = I(w - c) for a reference value c. The other forests describe the response where w is 0, and
    a column that never comes near 0, an age or a year, makes them an extrapolation and the fit poor,
    with nothing in the chain to show it. A dose or an exposure, where 0 means none, is left as it is.
    The forest enters as w (a f(x)), and sd is the size of a f(x), the change in the response for one
    unit of w; the term at an observation has size |w| sd. The column is used as given: it is never
    centred and never divided, so a stated sd is per unit of w, and the same sd with w in grams and in
    milligrams states two priors a thousand times apart. With several columns each has its own
    coefficient over the one forest, and sd is one number for all of them.

    Where the coefficient is held (amplitude = fixed()) "size" is the prior standard deviation, exactly.
    Where it is drawn the size is itself uncertain, and sd is its prior median.

    Not stated, sd is 2 sd(y) for a forest with no basis under a continuous response and one unit of the
    link's error otherwise; sqrt(2 / K) of sd(y), or of that unit, for a factor in a model of K forests;
    and for numbers that same size for one standard deviation of each column,
    sqrt(2 / K) sd(y) / sd(w) per unit. A default fit is then the same fit whatever units w is in.
    sd(w) is taken once, from the rows the sampler is created with, rows of weight zero included, and
    kept with the data: replacing the column later does not move the prior. A caller who replaces a
    column and wants the default of the new one states it,
        p <- sampler$getLeafPrior(2)
        sampler$setLeafPrior(forests = list(forest(), forest(sd = unname(p$sd * p$basis.scale / sd(w.new)))))
    A constant column has no standard deviation and needs sd stated; so does a basis written 1 + w.

    A treatment belongs in a factor or a logical vector, where sd is the treatment effect. Written as a
    0/1 number it is scaled as any number: sd(z) is 0.5, 0.4 and 0.22 at 50, 20 and 5 percent treated, so
    the default effect is 2, 2.5 and 4.6 times the one factor(z) gets;
    sd = sqrt(2 / K) * sd(y) on the number states the factor's. The indicator columns of a factor handed
    over as numbers, cbind(1 - z, z), are two coefficients: with sd stated the difference between the two
    has size 1.41 sd.

    A held coefficient cannot adapt to an sd that is too small. At the default it is sized for a slope of
    about one sd(y) per sd(w); where the effect may be larger, under "probit" and "logistic" above all,
    state sd or let the coefficient be drawn.

Until push 4 the item does not say that `print` shows what each forest resolved to; push 4 adds the sentence.

## Constraints

- A fit whose forests have no basis or a factor, coefficients drawn, draws what it draws on the base build,
  to the bit, in every push, stated sds included. In push 2 a fit with a drawn numeric multiplier does too.
- One scale for each column, from the record to the draw. No number per column is reduced to one number in
  R, in the forests' record, in the bridge or in the engine.
- The law is resolved in one function of the engine and nowhere in R: R states what the caller stated and
  the kind, and reports what the engine returns. The engine holds every default.
- A forest's law is resolved once for each construction, by the chain, and the combiner is handed the
  result; no second derivation can differ from the first.
- [`inst/include/dbarts/dbarts.h`](../../inst/include/dbarts/dbarts.h) is not edited and
  [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not move.
- The stored state is not edited: no block, no name, no encoding;
  [`stateFormatVersion`](../../src/R_interface_bartcore.cpp) stays.
  [`serializeGlue`](../../src/bartcore/combiner.hpp), [`restoreGlue`](../../src/bartcore/combiner.hpp) and
  [`glueIsValid`](../../src/bartcore/combiner.hpp) are not touched.
- `data@bases[[f]]` is not touched: no column is centred, divided or given an attribute.
- The unit, the conversion of a stated sd and what is recorded beside the anchor are forest-sd-unit's and
  are carried, not rebuilt. A default never passes through the unit on its way in.
- The kind is read from the one function forest-kind-by-class adds, and its swap table stands; this slice
  adds one row to it, the width of a default numeric forest.
- The lines ["m13"](../../benchmarks/R/mutation-battery.R) and m14 quote in
  [`drawForestAmplitude`](../../src/bartcore/combiner.hpp) stay as written; the entries that quote lines this
  slice rewrites (["m45"](../../benchmarks/R/mutation-battery.R) and what forest-sd-unit left of m46 and m47)
  move with them.
- Nothing of `updateBasisScale`, of `leaf.prior` on `forest()`, of a heavy-tailed coefficient or of an sd
  for each column rides along; each stays the refusal it is.
- Each push leaves the help saying what the code does. Base R calls stay within DESCRIPTION's R floor.

## Pushes

Four, each gated on its own and each a coherent tip.

1. Statements. The engine is handed what was stated (a kind, an sd or none, a hold) and holds the rule and
   the defaults itself; the rule is still the tip's. NEUTRAL. A tip on which every change of shape (the
   forests' record, the kind's route from the data to the engine, the tests/cpp fixtures) is proved bit for
   bit before any arithmetic moves.
2. Held coefficients. The held values go by kind, a held forest's sd is exact, the held shapes are the
   defined ones, and `amplitude.prior.variance` goes. The numeric law is untouched, so dec-A171's refusal of
   one held numeric column stays. bartCause's edit the same day.
3. The numeric multiplier. The column as given, the normal coefficient, the default for one standard
   deviation of each column with its record, one scale a column in the engine, the reader and `extract` by
   column, one held numeric column. dec-A171's refusal is lifted. The design note is completed here.
4. The printed block. R only.

No tip between them fits a model nobody asked for: push 1 moves nothing; after push 2 a held numeric
column is still refused, not held under half a law; push 3 lands the numeric law, its reader and its help
together. Push 3's tip lacks the printed block for as long as push 4 takes; its help does not promise one.

## Steps

"Fails today" is what the base build does where the test expects otherwise. New names are for the reader's
sake. Calls are in push 3 of written-surface's spelling; sds are in the response's units. Fixture, unless
said: 150 rows with x1, x2, a 0/1 z, its factor zf, a three-level g, dose (sd 0.02), age (mean 50, sd 17),
an offset column and weights of which ten are 0; a response whose standard deviation is neither 1 nor its
range.

### Push 1: statements

1.1 The engine. [`ForestSpec`](../../src/bartcore/combiner.hpp) states a kind (nothing, levels, numbers),
    the sd forest-sd-unit added (not a number where none is stated), the coefficient's variance where one
    is stated (push 2 removes it) and the hold; it loses its four derived numbers. One function (`forestLaw`)
    turns a statement into the forest's own sd, its divisor and its coefficients' prior, and one
    (`defaultForestSd`) holds the three defaults; in this push both give the tip's numbers, the row norm and
    the starts by position included, with every expression of
    [`mapLeafScale`](../../src/bartcore/chain.hpp) in its order of operations. The chain resolves each
    forest once and hands [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp) the result.
    [`expandForestSpecs`](../../src/bartcore/combiner.hpp) states BCF's two forests the same way. Thrown: a
    forest with no basis handed one; a levels or numeric forest handed none. Every struct here is read by
    objects that do not track headers: `--preclean`.
    Tests, tests/cpp: the fixtures of [`testForestMapWriters`](../../tests/cpp/test_sampler.cpp),
    [`testBCFCalibrationMap`](../../tests/cpp/test_sampler.cpp),
    [`testForestCalibration`](../../tests/cpp/test_sampler.cpp) and the 100 other lines that set the four
    numbers are restated, and every comparison they make holds unchanged. New (`testForestLawStatements`):
    for K = 2 and 3, each kind, held and drawn, stated and not, under each family, the forest's leaf scale,
    its coefficients' variance and its half-Cauchy median equal the tip's expression written out with
    literals, to the bit; the two throws.
1.2 The bridge and the record. [`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp) reads four numbers
    a forest (tree count, base, power, the hold) and the optional variance;
    [`parseData`](../../src/R_interface_bartcore.cpp) gives each forest its kind from the data object: no
    basis, a basis with the levels forest-kind-by-class records, a basis without.
    [`forestParams`](../../R/model.R) writes the four and stops deriving; `defaultAmplitudePriorScale` goes,
    and its one checked cite, in bcf-latent-evidence.md, is marked `retired:`. A record of eight numbers a
    forest is a sampler saved before: refused with the [1] text, where forest-sd-unit read it as it was.
    [`bartcore_getLeafPrior`](../../src/R_interface_bartcore.cpp) reports the kind the engine holds.
    Tests, a new file test-forest-law-statements.R:
    - The kind, through the engine. Each of eight objects (a factor, an ordered factor, a character and a
      logical vector, a 0/1 number, dose, two columns, the indicator columns of a factor) at each of the
      seven creation doors of forest-kind-by-class: the kind the bridge reports is the kind that slice's
      function gives for the data object (a kind dropped between the data and the engine reads "numeric"
      for a factor at that door).
    - Saved before. A sampler whose record is rewritten by hand to eight numbers a forest is refused with
      the text at `copy()`, at a run after `saveRDS` and `readRDS`, and at
      `new("dbartsSampler", control, model, data)`; a `bart` fit that kept its sampler, at `predict`.
1.3 Repair. 72 assertions pin the stored numbers by position on the stand-in (test-bcf-family.R 46,
    test-forest-arguments.R 16, test-bcf-creation.R 8, test-multiforest-leaf-prior-writer.R 2), and
    forest-sd-unit adds some: each is restated against the four numbers or against the reader.
    forest-sd-unit's test of a record with no `sd` and no `unit` becomes the refusal above. Run the suite
    first and repair what it shows.
1.4 Records. docs/architecture.md where it lists the forests' record; the mutation battery's entries.
1.5 Mutations (Verification).

### Push 2: held coefficients

2.1 The engine. The held rows of "The rule": a held forest with no basis has own sd s, by the leaf scale
    and not the half-Cauchy median; a held two-level block has no divisor. A held block is set to its
    kind's value when the combiner is built, at any position, after
    [`rebuildAmplitudeLayout`](../../src/bartcore/combiner.hpp), which keeps the drawn starts; a held shape
    with no value throws there. In this push one numeric column has none.
    [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp) on a held forest with no basis restates the
    leaf scale. The coefficient's variance leaves [`ForestSpec`](../../src/bartcore/combiner.hpp) and
    [`AmplitudeSpec`](../../src/bartcore/combiner.hpp): a level's is 1/2.
    Tests, tests/cpp (`testHeldCoefficients`): three forests (none, a two-level block, a three-level block
    that is drawn) in each of their six orders, the first two held: the held values by kind at every position
    (fails today in four of the six); each held forest's leaf scale equal to its sd over the unit times the
    anchor, to 4 ulp, with sds that differ from each other, from the unit and from 1 (fails today: the
    divisor, and a forest with no basis at the anchor); the throws for a held three-level block and a held
    numeric block of one and of two columns; a swap of a held two-level block to another of its width leaves
    (0, 1); a state stored from a twin that draws installs and leaves (0, 1); a write to a held forest with no
    basis gives the leaf scale of one created at that sd.
2.2 R. [`forest`](../../R/model.R) loses `amplitude.prior.variance`, and
    [`validateForestKnobs`](../../R/model.R) and [`resolveForests`](../../R/model.R) their lines for it. The
    held shapes are refused by name where the bases are in hand after `subset`
    ([`resolveSamplerSpec`](../../R/spec.R)), in the [2] texts, the forest named by its label;
    [`validateForestAmplitude`](../../R/model.R) takes the reworded `fixed(2)` text. The reader
    ([`reportLeafPrior`](../../R/dbarts.R)) gains `amplitude` (`fixed()` or absent), takes the four
    `prior.sd.of` values of a forest with no basis or a factor, and loses `amplitude.prior.variance`.
    Tests, a new file test-forest-held.R:
    - By kind, at every door and position. `amplitude = fixed()` on a forest with no basis, on `factor(z)`
      and on a logical, in a formula term, a `forests` list on a formula and on a matrix, `dbartsSpec()` and
      a data object's `bases`; as the first, second and third forest, and on a forest with no basis written
      second: `getForestAmplitudes()` is 1, or (0, 1), before and after 50 sweeps (fails today: (1, 1) third
      and first; 0 for a plain forest second).
    - Exact. The reader's `k.scale` of a held forest equals its sd: stated 0.3 and default on a forest with
      no basis, stated 0.6 and default on a two-level factor, under each family, to 1e-12 (fails today: the
      response's scale, and s / 0.674).
    - The size, from the engine's prior draws (`at_home`, 5000 sweeps a row): the second level's term of a
      held factor and a held forest with no basis have the stated sd within 4 percent.
    - Refused: a held factor of three levels, a held `cbind(1 - z, z)` and a held `dose + age`, each with
      its text at three doors (fails today: held at (0, 1, 1) and (0, 1)); `amplitude.prior.variance` is
      R's unused-argument error at both doors and through `do.call`.
    - Not moved: a model with one numeric column drawn beside a held factor keeps the base build's row norm
      (`k.scale` of the numeric forest to the bit), and dec-A171's refusal stands at creation.
2.3 The gates.
    - A new script, benchmarks/R/held-coefficient-exact.R: one binary predictor, one tree a forest, sigma
      held, no quadrature. Arms: a forest with no basis held at sds of 0.3 and 2 of the unit; a two-level
      factor held at 0.6, at 0.2 and at its default; both held. Each compares the posterior mean of each
      forest's term in each cell and the probability that each tree splits with the closed form written
      from "The rule"; in `quick` 8 seeds of 100000 kept sweeps, tolerance 0.01. Added to
      [exact-gates.yaml](../../.github/workflows/exact-gates.yaml)'s list.
    - The four BCF gates: each oracle takes the treatment forest's sd with no 0.674 in an arm that holds
      its coefficients and the prognostic forest's sd as stated in an arm that holds a, and the
      `amplitude.prior.variance` argument goes. [bcf-latent-exact.R](../../benchmarks/R/bcf-latent-exact.R)
      is the one that fails with its oracle left as it is (Context): run it once so on the push's build,
      six arms failing, and then with the oracle restated, none.
    - [sbc.R](../../benchmarks/R/sbc.R): the argument goes, and its fixed-glue arm draws its truth at the
      held scales; the arm constructs and runs its smallest setting.
2.4 The baseline. Run the BCF compare against `bcf-equivalence-1b7d730c.rds`: 14 scenarios identical, with
    no `max |z|` line, and ["glue_toggle"](../../benchmarks/R/bcf-equivalence.R) not. Re-record; the
    MANIFEST row names the oracle (rule P17), which is an identity and a gate: the re-recorded draws equal,
    to 1e-10, those the base build gives the same scenario with its treatment forest stating 0.674 of the
    default (pair row 33's identity, run for the scenario itself); and the held-factor arms of
    held-coefficient-exact.R and of bcf-latent-exact.R pass, with their gaps.
2.5 Respell. `amplitude.prior.variance` on 31 lines of 6 test files: where it states 0.5 the argument is
    dropped and the fit is the same; where it states another value the test pinned the tip's law and goes
    with it, one pin of R's error taking its place. Five benchmark scripts drop the argument. The pin of
    ["would hold the forest at zero"](../../inst/tinytest/test-forest-arguments.R) takes the shortened text.
2.6 Help and records. man/forest.Rd: the usage, the `amplitude.prior.variance` item removed, the
    `amplitude` item (the three held values, two in this push; that a held coefficient makes `sd` exact),
    the `sd` item's sentences for a held forest, the Details paragraph on the budget.
    [`dbartsSampler$getLeafPrior`](../../man/dbartsSampler-class.Rd) with its docstring. A new
    docs/design/forest-multiplier-law.md with its index row: the rule's held half, the changed sequences
    with their oracles, the re-recorded scenario. docs/design/bcf.md and multiplier-combiner.md where they
    give the 0.674 to a held forest or name the variance. TODO: `forest-prior-args`.
2.7 bartCause, same day (its own commit on dbarts-1.0; the kind-by-class slice's line,
    `basis <- factor(as.integer(z), levels = 0:1)`, must already be in, or `update.b = FALSE` is refused as
    two held numeric columns). In R/bcf.R: `b.prior.variance = 0.5` leaves the formals of `fitBCF` and of
    `bcf`; `amplitude.prior.variance = b.prior.variance` leaves the treatment forest's call;
    `b.prior.variance = b.prior.variance` leaves the call of `fitBCF`; and beside its refusal of `forests`
    and `bases`, `bcf` refuses the argument by name ("'b.prior.variance' is no longer an argument: the
    treatment coefficients' prior variance is fixed; state the effect's size with 'sd.moderate'"), since
    its dots would otherwise take it in silence. man/bcf.Rd: the usage line and the item go; `update.a,
    update.b` reads "when FALSE the coefficient is held, at 1 for the prognostic forest and at 0 and 1 for
    control and treated, and sd.control or sd.moderate is then that forest's prior standard deviation,
    exactly". tests/testthat/test-14-bcf.R and test-03-responseFit.R: the two hand-built samplers drop
    `amplitude.prior.variance = 0.5`, the same fit. A `bcf()` fit with `update.a = FALSE` or
    `update.b = FALSE` moves; no test there pins a draw under either.
2.8 Mutations (Verification).

### Push 3: the numeric multiplier

3.1 The scale. One function of the engine (`basisColumnScales`) gives each column's sample standard
    deviation over every row, n - 1, the mean taken in two passes; a column with none (constant, or fewer
    than two rows) throws. Tests, tests/cpp: columns of sds 0.25 and 8 give exactly those; a 0/1 column at
    20 percent; a column with half its rows zero equals the literal over ALL rows; 1e8 + 1e-8 k within
    1e-6 of the two-pass literal (fails with one pass: 18 percent); a constant column and one row throw.
3.2 The law. [`ForestSpec`](../../src/bartcore/combiner.hpp) gains a stated sd for each column (`columnSds`,
    empty for one number on all) and a scale for each column (`basisScale`, empty for "take it").
    `forestLaw`'s numeric rows are "The rule"'s: s_j stated, or the default over the scale handed or taken;
    the forest's own sd from s_1; each coefficient's prior sd s_j / s_1; one column held at 1. A ratio whose
    square is not representable throws. [`ForestAmplitudePrior`](../../src/bartcore/combiner.hpp) gains the
    scale of each coordinate, and [`drawForestAmplitude`](../../src/bartcore/combiner.hpp) divides each
    coordinate's prior precision by its square on a line of its own; with none the block draws as today. A
    drawn numeric block starts at its scales. [`basisRowNorm`](../../src/bartcore/chain.hpp) goes with
    everything that read it. One function of the chain (`applyForestLaw`) sets a forest's leaf scale and
    the combiner's scales from the sds in force, and the constructor and
    [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp) both call it; the write marks the forest stated
    and, where the number is the one in force on every column, moves nothing.
    [`Chain::setForestBasis`](../../src/bartcore/chain.hpp) installs the basis and derives nothing, and
    refuses another width on a numeric forest that is not stated.
    [`ForestCalibration`](../../src/bartcore/chain.hpp) reports, for each column, the sd in force as stated
    or resolved, the scale held, whether the sd was stated, and the scale the COMBINER holds for the draw.
    Tests, tests/cpp (`testNumericMultiplierLaw`):
    - Literals. K = 3 under each family, drawn and held: a forest with no basis, a two-level block and a
      numeric forest; then two numeric columns of sds 0.25 and 8 with no sd stated. Each forest's own sd,
      its divisor, each coefficient's prior sd and start, to 4 ulp (the powers of two make the defaults
      exact).
    - One scale a column (dec-B272). The two-column forest once with no sd and its scales taken, once with
      the two sds that default resolves to stated as `columnSds`: identical draws over 50 sweeps. With every
      row masked, 5000 sweeps: column j's a_j F has the quartiles of a one-column forest stated at s_j,
      within 3 percent. With the forests' structure frozen, each coefficient's conditional mean and
      variance over 4000 draws against the closed form whose prior precision is 1 / (s_j / s_1)^2.
    - The twin. w against 1024 w, and (w, v) against (1024 w, v / 64), default: identical draws, drawn and
      held, gaussian and probit.
    - One home. After construction, after a write of one number on the two-column default forest, and
      after a swap: the combiner's scale of each column times the forest's own sd times its divisor equals
      the sd the reader reports, to 4 ulp; and after the write 20 sweeps are identical to a twin constructed
      stated at that number and given the same state.
    - Paths in the engine. A chain made again with the first one's scales handed back has its leaf scales
      and scales to the bit, on a basis since swapped to ten times the column; made again with none handed,
      it differs. The throws.
3.3 The bridge. [`parseData`](../../src/R_interface_bartcore.cpp) reads `data@basis.scale` through the
    guard an absent slot needs; [`applyForestBases`](../../src/R_interface_bartcore.cpp) hands a numeric
    forest that states no sd its entry, checked for width and sign there and only there.
    [`bartcore_getLeafPrior`](../../src/R_interface_bartcore.cpp) carries the four per-column vectors and
    the stated flag; [`bartcore_setForestBasis`](../../src/R_interface_bartcore.cpp) raises the engine's
    refusal.
3.4 The record. A slot `basis.scale` on [`dbartsData`](../../R/A_class.R), beside the levels
    forest-kind-by-class adds: `NULL`, or a list with one entry a forest, `NULL` or one number a column;
    read through an accessor that answers `NULL` for an object saved without the slot
    ([`dataRowNames`](../../R/data.R)'s pattern). Written in the sampler's `initialize`, where the anchor is
    recorded, from what the engine reports, for a forest whose entry the object lacks; dropped for a forest
    whose basis a declaration replaces on a data object ([`resolveSamplerSpec`](../../R/spec.R)), and for a
    forest whose swap changes its width ([`setForestBasis`](../../R/dbarts.R)).
    Tests, a new file test-forest-scale-record.R. Each path asserts the recorded scale, the reader's `sd`
    and the forest's `k.scale * 0.674` against literals written in the test, with a dose whose sd is not 1:
    - Which rows. With `subset`, a dropped missing response, weights with zeros and a mask set after
      creation: the scale is `sd(dose[kept])` to 1e-12, `kept` written in the test and
      including the rows of weight zero; and a default forest's `sd` is `sqrt(2 / K) * U / sd(dose[kept])`.
    - Once. After `$setForestBasis(2, dose / 10)`: the record, `sd`, `k.scale` and the coefficients are
      `identical()` to what they were (fails today: the sd follows the column). Then `copy()`; `saveRDS`,
      `readRDS` and a first use; `new("dbartsSampler", control, model, data)`: each `identical()` in record,
      `sd` and `k.scale` to the sampler (fails today: ten times off). With the record cleared by hand, `new()`
      takes the scale of the column in force.
    - Other data. `dbartsSpec()` over `sampler$data` with a declaration that gives the forest
      `100 * dose`: the record and the sd are the new column's (a kept record gives 100 times the sd). A
      `forests` list and a control taken from a sampler, given to a fit on other data: its own scale. Two
      fits from one `dbartsData()` object the caller holds: equal records, and the caller's object carries
      none.
    - Carried and wrong. A record of another width or with a zero, on a forest that states no sd: refused
      with the bridge's two texts; on a forest that states one: created.
    - Stated, then wider. A default forest written `sd = 0.5`, swapped to two columns, copied, reloaded:
      `sd` is `c(0.5, 0.5)` at each, the record's entry is gone and both run (a stale record stops the
      copy).
    - Not a path. `predict` at new rows equals the coefficients times the new basis times the forests' fits,
      by `type = "forest"`, to 1e-12, with a new dose a thousand times as wide; a state from a sampler
      created on another column installs and the recipient's record and `sd` are its own; `setResponse`,
      `setOffset`, `setWeights`, `setPredictor` and `setActiveRows` leave the record and `sd` `identical()`.
3.5 Refusals and what is accepted. Where the bases are in hand after `subset`: a constant column on a
    forest that states no sd, and columns too far apart, by their [3] texts. In
    [`setForestBasis`](../../R/dbarts.R), in the kind-by-class slice's table: another width for a numeric
    forest at its default. [`samplePriorPredictive`](../../R/dbarts.R) refuses several forests at its door.
    dec-A171's refusal goes from [`resolveSamplerSpec`](../../R/spec.R), and `refuseHeldOneColumn` with it:
    one numeric column is held at 1. Tests: a constant column, and a column constant on the rows `subset`
    keeps, refused at four doors and created with `sd = 2`; `basis = 1 + dose` refused and created with an
    sd; the width refusal, with the sampler, its record and the next 5 sweeps `identical()` to a twin's
    afterwards, and the same swap accepted after a write; the six tests of dec-A171's refusal turned over:
    `amplitude = fixed()` on dose is created at every door, held at 1 as the second, third and first
    forest, with `k.scale` equal to the sd (stated, and `sqrt(2 / K) * U / sd(dose)`); swapped to another
    column it is still held at 1 with the same `k.scale`, and swapped to two columns it is refused by the
    kind-by-class slice's text; a held 0/1 number and the held factor of the same column, at one stated sd,
    have `identical()` train draws.
3.6 The reader, the writer and `extract`. [`reportLeafPrior`](../../R/dbarts.R) builds the entries of "The
    reader, `extract` and the printed block"; [`resolveForestSpreads`](../../R/dbarts.R) and
    [`writeForestSpreads`](../../R/dbarts.R) keep the writer's form, one number a forest;
    [`extractParameter`](../../R/generics.R) returns the list, named by the fit's labels, and refuses a fit
    whose stored entries have no `sd` with the [3] text;
    [`leafPriorIsDrawn`](../../R/bart.R) and [`defaultLeafScaleVars`](../../R/diagnostics.R) read
    `multiplier` where they read `basis.row.norm`. Tests, a new file test-forest-leaf-prior-reader.R:
    - The literals of the table, for a `bart` fit of four `forest()` terms and for the sampler, `sd.stated`
      before and after a write.
    - Through the engine. For every forest of that model and of its held twin, after creation, a write, a
      swap, a copy and a reload: the bridge's scale of each column as the combiner holds it, times
      `k.scale`, times 0.674 where the coefficient is drawn, equals `sd` to 1e-12; for a drawn forest with
      no basis the half-Cauchy median times the unit does.
    - Read then write. `s$setLeafPrior(forests = lapply(s$getLeafPrior(), function(p) p$leaf.prior))` on a
      default and on a stated sampler leaves `sd`, `sd.stated` and the next 20 sweeps `identical()`. A
      default numeric forest written `forest(sd = unname(p$sd))` has `sd.stated` `TRUE` and the next 20
      sweeps `identical()` to an untouched twin's.
    - `extract`: `identical()` to the list of the kept sampler's `sd` entries, named by label; by a label
      and by a position a list of one; a fit with no kept sampler the same from what it stored; a stored
      list with its `sd` entries removed by hand refused with the text.
3.7 The gates.
    - benchmarks/R/numeric-multiplier-exact.R, the eleven arms of Context, the oracle written from "The
      rule" with R's `sd()` over every row the sampler holds and nothing read off the sampler but the unit
      and the shift. Functionals: the posterior mean of the term in each cell at each value of the basis;
      the probability that the multiplied tree splits; for one column the probability that the coefficient
      is within 0.5, 1 and 2; for two the second moments of the standardized coefficients, on a grid placed
      by two earlier passes. In `quick` 8 seeds of 100000 kept sweeps, tolerance 0.01; in full 16 of 300000,
      tolerance 0.005. Added to the workflow's list.
    - The prior from the engine, tinytest, `at_home` (5000 sweeps a row under a flat likelihood): each
      numeric row of "The rule" within 6 percent, default and stated, one column and two, a column times
      1000, an age, a 0/1 number at 20 percent, held, under each family.
    - The units twin, tinytest, always on: dose against 1024 dose, and two columns each by its own power of
      two, default, drawn and held, gaussian and probit: identical draws.
3.8 Respell and repair. Run the suite first. Expected from the stand-in: the 13 assertions that pin the
    row norm and a coefficient's variance (["medianRowNorm"](../../inst/tinytest/test-bcf-family.R), 12 of
    them in that file) are rewritten as pins of "The rule" or go with it; its 27 creations on a constant
    column, and the one each in test-formula-terms.R and test-sampler-residuals.R, state an sd; the 28 pins
    of the reader's entries and labels in 3 files
    (["forest total"](../../inst/tinytest/test-calibration-midchain.R)) take the new ones; the assertion
    forest-kind-by-class leaves, that a factor's indicator columns handed over as numbers draw what the
    factor draws, is turned around, with the 1.41. Each test that only needs some numeric forest is left,
    and now fits the law.
3.9 Help and records. man/forest.Rd: the `sd` item above, `amplitude` with one numeric column, `basis`
    where it speaks of a constant column. man/dbartsData.Rd: the slot.
    [`dbartsSampler$setForestBasis`](../../man/dbartsSampler-class.Rd),
    [`dbartsSampler$getLeafPrior`](../../man/dbartsSampler-class.Rd) and
    [`dbartsSampler$setLeafPrior`](../../man/dbartsSampler-class.Rd) with their docstrings: a swap moves no
    prior; the entries; a write states every column. The `extract` page's `leaf.prior.sd` item.
    docs/design/forest-multiplier-law.md completed: the rule, the scale and its record on every path, what
    a state and a saved sampler hold, the changed sequences with their oracles, and the four additions
    left with the reason each is one. docs/design/multiplier-combiner.md rewritten around the rule;
    nameable-calibration.md, public-surface.md, state-not-model.md and docs/architecture.md where they
    describe the row norm or the data object's slots. TODO: `forest-prior-args`, and the entries of "Out of
    scope".
3.10 Mutations (Verification).

### Push 4: the printed block

4.1 [`fitDescriptors`](../../R/bart.R) stores, for each numeric basis column, its mean and standard
    deviation as the run ended, so a fit prints with no sampler.
4.2 One function (`forestPriorLines`) builds the block from the reader's entries and those two numbers;
    [`fitSynopsis`](../../R/generics.R) prints it after the tree counts forest-defaults-by-kind prints, `show`
    of a sampler prints it, and a sampler created with `verbose = TRUE` prints it once. A fit of one forest
    and a fit stored without the entries print as they do today.
4.3 Tests, a new file test-forest-print.R: the block of "The reader, `extract` and the printed block" by
    exact string for the model of four forests, at `print`, `show` and creation; a stated forest, a held
    one, probit and logistic headers; after a swap to dose / 10 the column's sd is the new one and the
    line ends with the recorded one; a column centred by `scale()` prints mean 0; every number in the block
    read back from the text equals the reader's, or R's `mean` and `sd` of `data@bases`, to the digits
    printed; nothing is printed for one forest; no warning anywhere, for an age and for a 0/1 number.
4.4 Help: the `sd` item's sentence that `print` shows what each forest resolved to, with the two lines to
    look for; man/bart.Rd where it describes `print`.
4.5 Mutations (Verification).

## Verification

Every push, against the slice's own library (`R CMD INSTALL --preclean -l <lib> .`, `R_LIBS=<lib>` on every
call; check `dbarts:::buildInfo()$mode` and that the install postdates the source), run in series, at most
two cores (`MAKEFLAGS=-j2`, `EQUIVALENCE_CORES=2`):

- `cd tests/cpp && make && ./test_bartcore`, and again built with `OPT="-O2 -g -fsanitize=address,undefined"`
  under `ASAN_OPTIONS=detect_container_overflow=0`; the R-loaded path under the address sanitizer for the
  push's new test files ([Gate hygiene](README.md#gate-hygiene) gives the commands). Push 4 touches nothing
  under src/: tests/cpp unchanged and passing.
- The full tinytest suite on the shipped build, in one process, counted file by file: no failure, no file
  stopping, and at least the base build's count plus the new files' assertions less those removed with a
  law; the landing note gives the figures.
- On a reference build (`--preclean --configure-args=--enable-reference-build`): the four
  `test-reproducibility-*.R` files pass unchanged, and the three compares are bitwise, every scenario
  reporting identical draws, counted scenario by scenario with no `max |z|` line: 55 against
  `equivalence-1b7d730c.rds`, 15 against the BCF baseline in force, 11 against
  `multinomial-equivalence-80b1c8d4.rds`. In push 2 the BCF compare is first run against the old baseline,
  14 identical and `glue_toggle` not, and then against the re-recorded one, 15. `glue_toggle` fails the
  statistical mode against the old baseline, as a changed posterior must (111 of 209 summaries beyond
  z = 3 on the stand-in); what shows its new values right is step 2.4's identity.
- Every gate [exact-gates.yaml](../../.github/workflows/exact-gates.yaml) lists, in `quick`, each on its own
  exit status; in pushes 2 and 3 the four BCF scripts and the two new ones also in full.
- `inst/include/dbarts/dbarts.h` has no diff and `tools/check-api-hash.sh` passes.
- The pair script (below), old side on the base build, new side on the push's.
- `Rscript -e 'lintr::lint_package()'`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status;
  `pkgdown::check_pkgdown(".")`; `R CMD build` with every vignette rebuilt and `R CMD check --as-cran` on a
  tarball from a clean copy.
- The consumers, each suite whole against a private install of the push, none failing: bartCause on
  dbarts-1.0 (1412 expectations at its last run; with step 2.7's edit from push 2, and once without it on
  push 2, failing only where `amplitude.prior.variance` is passed), stan4bart on bartcore (582), treatSens
  on dbarts-1.0 (306), bairrtt on main (207), the last three unedited.
- A hot path is touched in push 3, one division a coordinate where a block's prior is seeded:
  `bench-sampler.R compare` on a quiet machine is the maintainer's to run.

The pair script. Each row is one seeded model on both builds: sampler fits compare the train draws, sigma
and the coefficients, `bart` fits `yhat.train` and sigma. Run for this plan on the tip against the stand-in,
which is the four pushes at once: rows 01 to 20 and 31 identical, rows 21 to 30 and 32 differing, row 35
created.

| | model | push 1 | push 2 | push 3 |
|---|---|---|---|---|
| 01 to 04 | two forests, default; an sd on both; three with a factor third; 8 levels with an sd | identical | identical | identical |
| 05 to 10 | probit; logistic; weights and an offset; a term under `subset`; `bart` with two chains; a logical | identical | identical | identical |
| 11 to 15 | a swap to another factor; a write; a copy; a reload; `setResponse` | identical | identical | identical |
| 16 to 20 | no forest without a basis; forest weights and a mask; four forests; one forest; counts and a tree prior on a factor forest | identical | identical | identical |
| 21 to 26 | a numeric column, default and stated; a 0/1 number; two columns; indicator columns as numbers; a factor beside a number | identical | identical | differ |
| 27 to 30 | a held factor; a held factor third; a held forest with no basis, default and stated | identical | differ | as push 2 |
| 31 | a held forest with no basis, default, probit | identical | identical | identical |
| 32 | a factor with `amplitude.prior.variance = 2` | identical | refused | refused |
| 33, 34 | against the base build's own spelling of the same prior: a held two-level factor, second forest, at s, the base build stating 0.674 s; a held forest with no basis at one unit, the base build stating any sd | | equal to 1e-10 | |
| 35 | a held numeric column | refused | refused | created |

Rows 33 and 34 are two changed sequences against a prior the base build can state: run for this plan on
the tip and the stand-in, four such pairs under gaussian, probit and logistic were identical. Where the base
build has no such spelling (rows 21 to 26, and a held factor outside the second position) the oracle is the
exact gate and the tests/cpp literals.

Mutations, each expected to fail the named test and no gate before it. Apply, install with `--preclean`
where it is the engine's, run, record the failing count, revert, `touch` the file.

- push 1, the default of a forest with a basis computed for K - 1 forests: step 1.1's literals at K = 3;
- push 1, a forest with no basis given the multiplied forest's default: step 1.1, and pair row 01;
- push 1, the kind not carried from the data door (every basis numeric there): step 1.2's kind test;
- push 1, the hold flag read from the wrong one of the four numbers: pair row 27, and the suite's held pins;
- push 1, the stated variance dropped: pair row 32;
- push 1, a record of eight numbers read as four: step 1.2's "saved before";
- push 2, a held value left by position: step 2.1's six orders, and step 2.2's third forest;
- push 2, a held forest with no basis second held at 0: step 2.2's plain forest written second;
- push 2, the 0.674 left on a held factor: step 2.1's leaf scale, step 2.2's `k.scale`, and
  held-coefficient-exact.R (0.10 and 0.21 at its two stated arms on the draft);
- push 2, a held forest with no basis left at the response's scale: the same three (0.033);
- push 2, its default taken as 1 under gaussian: step 2.2's default;
- push 2, a held factor of three levels held at (0, 1, 1): step 2.2's refusal, and step 2.1's throw;
- push 2, `amplitude.prior.variance` taken and ignored: step 2.2's unused-argument pin;
- push 2, a write to a held forest with no basis sent to the half-Cauchy median: step 2.1's write;
- push 2, a held block installed from a state: step 2.1's state;
- push 3, the default per unit (the scale read as 1): step 3.4's `sqrt(2 / K) * U / sd(dose[kept])`, and the
  gate (0.06 to 3.2);
- push 3, n for n - 1: step 3.1's 0.25 and 8, and step 3.4's 1e-12 (the gate cannot see it);
- push 3, the scale over the nonzero rows; over the rows of positive weight; before `subset`; with weights:
  step 3.1's half-zero column, step 3.4's "which rows", and the gate's weights and two-column arms;
- push 3, one pass for the mean: step 3.1's 1e8 column;
- push 3, the second column on the first one's scale, in the record's writer, in the bridge's read and in
  `forestLaw`, three mutations: step 3.2's literals, step 3.6's table, and the gate (1.36);
- push 3, the two scales exchanged; the precision divided by the scale and not its square: step 3.2's
  fixture and twin, and the gate (0.79);
- push 3, a numeric coefficient's variance left at 1/2; the row norm left in; the column centred: step
  3.2's literals, and the gate (0.11; 0.51);
- push 3, a stated sd divided by the scale: step 3.2's literals, and the gate's stated arms (0.14 to 0.65);
- push 3, the scale taken again at a re-creation; at a swap: step 3.4's "once";
- push 3, the record kept where a declaration replaces the basis: step 3.4's "other data";
- push 3, the record written back on every creation: the caller's object in step 3.4;
- push 3, an unread record checked: step 3.4's "carried and wrong" and "stated, then wider";
- push 3, a write that leaves the combiner's scales: step 3.2's "one home", and step 3.6's "through the
  engine";
- push 3, a write of the number in force skipped whole, the forest left at its default: step 3.6's
  `sd.stated`, and step 3.5's swap after a write;
- push 3, a numeric block started at 1: step 3.2's two-column twin;
- push 3, a held numeric column at 0: step 3.5, and the gate's held arms (0.27);
- push 3, the reader's `sd` taken from the statement with the engine's scales unread: step 3.6;
- push 3, `extract` left a vector; named by position: step 3.6;
- push 3, the constant-column check made before `subset`; the width refusal dropped; the prior-predictive
  refusal dropped: step 3.5;
- push 4, a column's printed sd read from the record; the per-sd figure from the column in force where the
  default came from another; a label by position; a mean not passed through `zapsmall`: step 4.3.

## NEWS

No new item: forests, `forest()`, its `sd` and `amplitude`, `$setForestBasis` and every reader here are new
in 1.0-0 and nothing released changes. The "Multi-forest models" item is reread and changed only where it
names `amplitude.prior.variance` or the 0.674.

## What this leaves for after the merge, and why each is an addition

Each requirement with what holds it; "test" is the step above that fails if it is broken.

| later | what it needs of 1.0 | held by |
|---|---|---|
| `$setForestBasis(updateBasisScale = TRUE)` (dec-B265, dec-B276) | the argument does not exist, so no call states it | R's unused-argument error, pinned in step 3.5 |
| | with no argument a swap keeps the scale: `FALSE` is today's meaning | step 3.4's "once" |
| | a record for each column that a call can rewrite, and the arithmetic to restate a forest from scales | the slot of step 3.4; `basisColumnScales` and `applyForestLaw` of 3.1 and 3.2 |
| | no state, no shipped signature: the update is a virtual of the engine's own facade and one more bridge argument | Constraints |
| `leaf.prior = normal(sd = )` on `forest()`, and `normal()` for the default back (dec-B266, dec-B276) | `leaf.prior` on `forest()` is R's unused-argument error | written-surface's pin |
| | the reader's `leaf.prior` already says "stated" or "the default", and `sd.stated` with it | step 3.6's table |
| | a stated forest keeps its record while its width stands, so a default put back has the creation's scale | step 3.4's "stated, then wider" pins when it goes |
| a heavy-tailed coefficient, `amplitude = student(3)` (dec-B262, dec-B271) | `amplitude` takes `fixed()` and nothing else | written-surface's pins |
| | one row of `forestLaw`: a drawn numeric or factor block under the forest with no basis's law | step 3.2 puts the law in one function |
| | the state already carries one variance a forest and installs it where a size is half-Cauchy | Context; no edit |
| | the reader has an `amplitude` entry to say it in | step 2.2 |
| an sd for each column, in any of the three spellings (dec-B272, after 1.0) | the engine takes one stated sd a column | step 3.2's fixture |
| | one `sd` on several columns is that size for one unit of each, now and then | "The rule" 1; step 3.7's two-column stated arm |
| | the reader, `extract` and `print` are by column from the start | steps 3.6, 4.3 |
| | no stored state and no compiled interface holds a forest's sd | Constraints |
| | the forests' record holds one number a forest; one a column is read beside it, a record of one number meaning all | nothing to build now |

## What waits on what

- On push 3 of written-surface: a basis's column names, which name the reader's `sd` and the printed
  lines; the labels that name the reader's and `extract`'s lists and the texts; `basis = 1 + dose`, a
  constant column this slice refuses without an sd; every call here is in its spelling.
- On forest-defaults-by-kind: a forest with no basis may stand anywhere, so step 2.2 holds one written
  second; selection by label in the reader; its line of tree counts, which the block follows. Until push 2
  lands, `amplitude = fixed()` on a forest with no basis that is not the first holds it at 0 (Context), on
  every tip from that slice's landing: that slice should refuse it, or this push follow it closely.
- On forest-sd-unit: the unit, the engine's stated sd and "sd in force", the recorded unit beside the
  anchor. Push 1 replaces what that slice leaves of the four numbers and turns its reading of an old record
  into a refusal. Push 3 changes what its reader returns for a default forest: `forest()` in `leaf.prior`
  and the number in `sd`, where that plan returns `forest(sd = <the default>)`; and a write of the number in
  force, which that plan skips, here marks the forest stated and moves nothing.
- On forest-kind-by-class: the record of a factor's levels and the one function that gives a kind; its swap
  table, which gains a row; its push 2, which makes every block that stands for a factor a factor, so that
  what moves in push 3 is a number that was meant as one; bartCause's factor line, before push 2.
- To recheck once each has landed: the stand-in is rebuilt on the new tip and the suite rerun, for the
  counts of Context (132, 72, 29, 31); the pair script's 20 rows; that no baseline scenario has gained a
  number or a hold; that `samplePriorPredictive` still cannot predict a sampler of several forests once
  push 3 of written-surface rebuilds a basis at new rows (if it can, it returns draws with the
  coefficients at their starting values, and step 3.5's refusal moves to push 1); the names of the
  functions cited here.
- With [leaf-conversions.md](leaf-conversions.md) and
  [cross-family-state-install.md](cross-family-state-install.md): they edit the state's reader and writer,
  `Chain`'s install and `setData` paths and the leaf models; this slice edits the K-forest constructor, the
  three forest setters, the amplitude combiner's constructor and draw, and the bridge's creation path, and
  no line of the state. Serial for the shared files, in either order, each with `--preclean`.
- With the state-frame arc's data record: `basis.scale` is a fourth record of its kind, data derived once
  and then held, and follows its rule (an object that carries one: used). This slice adds the slot and its
  writer itself and waits for nothing; whichever lands second folds the two read-backs into one. Two things
  that arc must not do to it: write it after a state install, which never moves it, or take it off the
  object a copy hands to creation, which reads it there; step 3.4's copy after a swap is the guard.
- With the control-migration arc: its first slice moves the forests' record to the model, whole, as push 1
  leaves it. Serial; this slice first keeps the posterior-changing gates away from slices gated bit for bit.

## Out of scope, and where it goes

- `updateBasisScale`, the long form on `forest()`, the heavy-tailed law: after the merge (TODO
  `forest-after-merge`, `coefficient-law-opt-in`). An sd for each column: after 1.0 (dec-B272).
- `setData` on a sampler of several forests: refused, as today. The rule for the day it opens goes in the
  design note: a new object's record is used, what it lacks is taken.
- Holding or releasing a coefficient on a live sampler: no setter exists; `setModel`'s refusal of a model
  that differs in it waits for the forests' record to be on the model.
- To TODO as new entries: a drawn coefficient starts by its forest's position (1, and 0 for the second
  forest's first), a start by kind being a change of draws for a model of three forests;
  `samplePriorPredictive` for several forests, which needs the coefficients drawn from their prior; a state
  whose held block differs from the recipient's installs under the standing rule, with no message; a
  consumer that creates a sampler of several forests through the shipped header gets a scale taken at each
  creation, there being no R object to record it on (none does).

## Calls made in planning

- Four pushes, cut between the statements, the held coefficients, the numeric law and the printed block.
  One push is about 7800 lines in one review. The first cut puts every change of shape under a bitwise gate;
  about 100 lines of the tip's law are written to be replaced. The second lets the held law land on the
  gates that exist and gives the numeric law a hold that already goes by kind; its cost is that one held
  numeric column stays refused for one push more.
- The exact gate is two scripts, one a push, where the design names one. The gaussian BCF gates cannot see
  the held law at their designs (Context), so the held law gets arms that can.
- A two-column arm, with a quadrature over two coefficients, where the design had none: it is what sees a
  scale collapsed or exchanged between columns, the fault this slice is most open to.
- The gate reads nothing off the sampler but the unit and the shift, and takes the scale from R's `sd()`.
  The alternative, reading the scale back, would pass any wrong scale.
- The reader's `leaf.prior` is the statement, a default as `forest()`, and `sd` the numbers (dec-B269's
  words); a default of two columns has two numbers and `forest(sd = )` takes one.
- A write of the number in force marks the forest stated. Skipped whole, a caller told to state an sd
  before a wider swap could not do it with the number the reader gives.
- A factor's entry in `extract` is one unnamed number, the difference between two levels: a level's own
  coefficient has that sd over 1.41, and the number under each level's name would be false.
- `extract`'s result is a list whatever `forest` selects. The alternative, a bare vector for one forest,
  makes the shape depend on the length of an argument.
- The reader says "factor" where the kind-by-class slice's function says "levels": one is the word the
  user wrote, the other names the columns.
- Gone from the reader: the four entries of the tip's channels. A variance of 1 for a forest would be false
  where two columns differ, and the divisor is in `prior.sd.of`'s "(prior median)".
- The scale is unweighted and counts rows of weight zero, as the response's scale does: one estimator for
  both, and the one R's `sd()` gives for the column the caller sees.
- The mean is taken in two passes: on a column far from zero one pass is 18 percent off R's `sd()`, and the
  printed line would name a scale the user cannot reproduce.
- A stated forest keeps its record until its width changes, and loses it then. Kept, a default put back
  later has the creation's scale; dropped on a width change, no stale entry stops a copy.
- A record is dropped where a declaration replaces a forest's basis on a data object. A hand-edited object
  keeps its record, as the state-frame design rules for its three: to take a new one, clear the entry.
- A sampler saved before push 1 is refused, not read as it was. Reading it would keep the row norm, the
  0.674 on a held forest and the held values by position alive in the engine for objects only a development
  build could have written.
- A state holds nothing of the law, so none is refused and none is re-read. A state stored on a build that
  held a forest's coefficients at other values installs its trees under the recipient's, by the standing
  rule; a check there would refuse what
  ["fixed amplitudes"](../../inst/tinytest/test-state-not-model.R) pins as taken, a state whose held block
  is not the recipient's.
- Drawn starts stay where they are for a forest with no basis or a factor. A start by kind would be the
  cleaner rule and would move the draws of every model of three forests.
- A held numeric coefficient with no sd takes the drawn default, as the design has it: right at ordinary
  slopes, short under probit when the slope is two or three times larger, which the help says.
- A held forest with no basis at its default has sd 2 sd(y) under gaussian, the drawn default, where the
  tip has sd(y): one default for a kind, held or drawn.
- A constant column's refusal names `sd =` and the help names `1 + w`: push 3 of written-surface accepts
  that basis, and here it needs an sd.
- The verbose summary gains the block from R; the engine's own summary, which prints the first forest's
  model, is left to the TODO entry forest-sd-unit adds.
- Measured for this plan on stand-ins, not on the slice: the law's 28 rows, the 20 pairs, the suite's 132
  and the BCF compare come from the design's prototype ported to the tip, which holds the kind as an
  attribute, the sd in units of L and no four-number record; the gate's figures from a draft run against
  the prototype's own build, in its spelling. Read and not run: the full gaussian compare but for its one
  scenario of several forests, the snapshot files, sanitizers, the consumers' suites.
