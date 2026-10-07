# forest-multiplier-law: what a forest's coefficient is, held or drawn, and what a numeric multiplier's sd is the size of

Status: PLANNED (dec-B253 and dec-B246 as revised; dec-B260 to dec-B265; dec-B269 as dec-B280 and dec-B282
revise it; dec-B271, dec-B272, dec-B275, dec-B276, dec-B281, dec-B282; dec-A171). Follows
[forest-kind-by-class.md](forest-kind-by-class.md), [forest-sd-unit.md](forest-sd-unit.md),
[forest-defaults-by-kind.md](forest-defaults-by-kind.md) and push 3 of
[written-surface.md](written-surface.md), none of which has landed. Amended 2026-10-07 after its blind
critique (verdict: build after corrections) and after dec-B281 and dec-B282, under which a basis is one
numeric column or a factor of two levels: three pushes where there were four, and nothing for several
columns.

One plan, three pushes. The engine's change of law and the R surface that states it are not planned
apart: at every tip the help must say what the engine does and the reader must report the number the
engine uses, and a plan for the engine alone would leave a tip that does neither.

agent: pushes 1 and 2 (the law): opus implementer for the engine, the bridge and the R code, sonnet for
the respelled tests, the benchmark scripts and the help once the code is fixed, opus reviewer told to
refute. Push 3 (the printed block): sonnet implementer, opus reviewer. The reason for opus on 1 and 2:
every slip is a prior off by a factor and no message. The scale of a column can be taken from the wrong
rows or taken twice on any of nine paths, and a held value can go by position on one path and by kind on
another; a stand-in for the slice, built with care, had three such faults (Context).
rng: stated per push and per call sequence, as [RNG classes and their gates](README.md#rng-classes-and-their-gates)
defines the classes.
- Push 1: POSTERIOR-CHANGING for a model that holds a coefficient, in two sequences: a held two-level
  factor (the forest's own sd is s where it was s / 0.674); a held forest with no basis that states an
  sd, or is left at its default under a gaussian response (its sd is s where it was the response's scale
  whatever s; under probit and logistic the default is the same number and the fit keeps its bits).
  NEUTRAL for every model whose coefficients are all drawn, a numeric multiplier included. Accepted
  where refused, forest-defaults-by-kind's two interim refusals being lifted: a held forest with no basis
  that is the second forest, at 1; a held two-level factor that is not the second, at 0 and 1. Refused
  where accepted: `amplitude.prior.variance`.
- Push 2, its first commit, the bound under which a row leaves a forest's update: bit for bit for every
  recorded scenario, every snapshot and every seeded pin of the suite, and gated so, alone. In law it
  changes a sweep in which some row's multiplier is under 2^-26 in absolute value and not under 2^-26
  of the forest's largest, or the reverse; that is every sweep of a numeric column whose values are
  below about 1e-8, to which the forest was blind, and otherwise a sweep where a coefficient has come
  within nine orders of magnitude of zero or of another level's.
- Push 2, the law: POSTERIOR-CHANGING for every model with a numeric multiplier, drawn, at the default
  sd or a stated one. NEUTRAL for every model whose forests have no basis or a factor, held or drawn,
  stated or not. Accepted where refused: a held coefficient on a numeric column. Refused where accepted:
  a numeric basis that is constant, or constant to rounding, on a forest that states no sd.
- Push 3: no draw moves.
- Recorded scenarios. Push 1 moves one, `glue_toggle` of the BCF equivalence baseline, which holds the
  treatment forest's coefficients, and re-records it. Push 2 moves none and adds four to that baseline,
  so that the law it builds occurs in a recorded draw. No other scenario of the three baselines holds a
  coefficient or multiplies a number, and the four snapshot files fit no model of several forests
  (searched; the BCF compare was run on a stand-in: Context).
window: pre-release, after forest-kind-by-class's second push and before the control-migration arc moves
the forests' record to the model. Serial with any other work in
[`ForestSpec`](../../src/bartcore/combiner.hpp), [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp),
the K-forest constructor of [`Chain`](../../src/bartcore/chain.hpp),
[`Chain::setForestBasis`](../../src/bartcore/chain.hpp), [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp),
[`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp), [`parseData`](../../src/R_interface_bartcore.cpp),
[`bartcore_getLeafPrior`](../../src/R_interface_bartcore.cpp), [`forestParams`](../../R/model.R), the
multi-forest block of [`resolveSamplerSpec`](../../R/spec.R), [`reportLeafPrior`](../../R/dbarts.R) or the
sampler's `initialize`; so serial with [leaf-conversions.md](leaf-conversions.md), in either order (What
waits on what). bartCause's edit lands the day push 1 does.
budget: ~5400 lines changed, upper figure 8900. By push: 1 ~2050 (upper 3380), 2 ~2830 (upper 4670),
3 ~520 (upper 850). By layer, with the upper figure: engine 370 (610); bridge 230 (380); tests/cpp 600
(990); R 740 (1220); tinytest 1720 (2840); benchmarks 800 (1320); help 430 (710); design note,
architecture, TODO and the two indexes 510 (830). As first planned the slice was 7800 lines in four
pushes. Gone from it: the neutral first push, folded into forest-sd-unit (about 950, the kind's route
staying here); everything for several columns and for a width that changes (about 1300); the refusals
of development-build objects and a second gate script (about 300). Added: the relative bound with its
twins (about 150), four scenarios (80), an arm of three trees (60), the tests that run always (40).
Past estimates in this arc ran 1.5 to 2 times low; the upper figure is 1.65 times the estimate.

## Goal

Each kind of forest has one law, and it is the law the help states. A forest's `sd` is the size of what it
contributes for one step of its multiplier: exactly where its coefficient is held, as the prior median where
it is drawn. A numeric column is used as it was given, never centred and never divided; its coefficient is
normal; a stated sd is per unit of the column; and with no sd stated the default is sqrt(2 / K) of the
response's scale for one standard deviation of the column, that standard deviation taken once, when the
sampler is created, and kept with the data. A held coefficient has the value its forest's kind gives it,
wherever the forest stands: 1 with no basis, 0 and 1 for the two levels of a factor, 1 for a numeric
column. `amplitude.prior.variance` is not an argument. A default fit is the same fit whatever units a
column is in, down to units in which its values are of order 1e-12. The reader, `extract` and `print`
give a numeric forest's number under its column's name, and it is the number the engine draws under.

## Context

Measured at d02fe72c (written-surface pushes 1 and 2 landed; push 3, forest-defaults-by-kind,
forest-sd-unit and forest-kind-by-class not) on the shipped build, R 4.6.1, unless a line says otherwise;
R/ and the combiner are byte for byte the tip's since. The spellings are the tip's, a tilde on a basis
written as code. L is the response's scale as the engine holds it: sd(y - offset) under gaussian, 1 under
probit, pi / sqrt(3) under logistic; forest-sd-unit makes it the unit an sd is divided by and changes
nothing else below. K is the number of forests. "Size" is the prior median of the absolute multiplier
times the prior standard deviation of the forest's own fit, from the engine's prior draws under a flat
likelihood (4000 sweeps a row, Monte Carlo error about 0.02); a held multiplier is exact. Fixture: 300
rows; x1, x2; a 0/1 z; w with standard deviation 3.15; an age with mean 50. Rows for a factor of three
or eight levels and for two numeric columns were measured too; both shapes are refused before this slice
(dec-B281, dec-B282) and are left out.

- The suite. At the tip 232 files; on push 3's build (4eeaf03d, traced 2026-10-07) 233, of which 230
  without the three that ask for more than two threads give 17121 results, none failing.
- The law today. s is the stated or default number, c the median nonzero absolute value of the column.

  | forest | coefficient | the forest's own sd | what was measured, over s L |
  |---|---|---|---|
  | no basis, drawn | N(0, v), v inverse-gamma: its size half-Cauchy, median s (measured 2.05 at s = 2; q90 / q50 6.6, the law's 6.3) | L (0.999) | a F: 1.02 |
  | no basis, held | 1 | L, whatever s (0.995) | F: 1.42 at s = 0.7; 0.50 at the gaussian default 2; 0.99 at the latent default 1 |
  | two-level factor, drawn | N(0, 1/2) each (sd 0.71) | s L / 0.674 (1.48 at s = 1) | (a_2 - a_1) F: 0.99 to 1.01 |
  | two-level factor, held, second forest | (0, 1) | s L / 0.674 | (a_2 - a_1) F: 1.49 |
  | one numeric column, drawn | N(0, 1/2) | s L / (0.674 c) | a F per unit: 0.72 / c; per sd of the column 1.12 for w and 0.24 for the age |
  | one numeric column, held | | | refused (dec-A171) |

  The defaults: 2 for a forest with no basis under gaussian and 1 under probit and logistic; sqrt(2 / K) for
  one with a basis, a factor and a number alike. The same figures hold under probit and logistic (nine rows
  measured, each the gaussian row's to Monte Carlo error).
- `amplitude.prior.variance = v`, on a forest with a basis: the coefficients are N(0, v) and the forest's own
  sd does not move, so the size is sqrt(2 v) times the table's: 1.95 at v = 2 and 0.52 at v = 1/8 on a factor
  (2.04 with `sd = 0.7` beside it), 1.03 / c at v = 1 on a number. Nothing on a held forest. Refused on a
  forest with no basis. Stated on 31 lines of 6 test files, in 5 benchmark scripts and by bartCause.
- A held value goes by position. [`AmplitudeState`](../../src/bartcore/combiner.hpp) starts at 1 for the
  first forest and (0, 1) for the second, and [`rebuildAmplitudeLayout`](../../src/bartcore/combiner.hpp)
  carries those by position and fills every other coefficient with 1. Read off the engine, and again off
  push 3's build: a two-level factor held as the third forest, or as the first where no forest is
  without a basis, is (1, 1): no contrast, a function added on every row (size of the difference 0.000);
  a forest with no basis held as the SECOND forest is held at 0. A drawn coefficient starts at the same
  values. forest-defaults-by-kind refuses those two held shapes, by width and position, until this
  slice's first push; the suite holds 49 forests on push 3's build, every one a forest with no basis
  first or two columns second.
- The row norm. c is taken over the rows the sampler holds, after `subset` (1.976 for rows 1 to 200, R's
  median of those rows), over the rows where the basis is not zero, unweighted; rows of weight zero and rows
  masked later count. It is taken again at every swap: `$setForestBasis(2, 3 * w)` divides the forest's own sd
  by 3 (2.493 to 0.831) and leaves the coefficient. `copy()`, a reload and
  `new("dbartsSampler", control, model, data)` take it from the column then in force, so after a swap to
  w / 10 all three have ten times the creation's sd. Nothing of it is recorded.
- A row leaves a forest's update when its multiplier is under 2^-26 in absolute value
  ([`zeroMultiplierTolerance`](../../src/bartcore/combiner.hpp), in
  [`AmplitudeForestCombiner::formForestResponse`](../../src/bartcore/combiner.hpp)). The bound was written
  for BCF, whose multipliers are of order 1; a numeric multiplier is a coefficient of order 1 times the
  column as given, so the bound is in the column's units. Run by the critic and again by the previous
  planner, w against 2^-k w with one seed, on the tip and on the stand-in: identical draws down to
  2^-10 and 2^-16, not at 2^-20, and at 2^-30 and 2^-40 the forest sees no row (mean sigma 0.497
  becomes 2.43 on the tip, 0.490 becomes 2.52 on the stand-in; a held column on the stand-in 0.492
  becomes 2.34), with no message. On a copy of the stand-in with the bound taken relative to the
  forest's largest absolute multiplier in the sweep, six lines: identical to the bit at every power to
  2^-40, drawn and held, and six seeded fits of forests with no basis or a factor still `identical()` to
  the tip. Not run on that change by anyone: the three baseline compares, and the bench.
- Other paths. `setData` is refused on every sampler whose forests carry coefficients, and
  `setResponse(updateScale = TRUE)` too; `setResponse(updateScale = FALSE)`, `setOffset`, `setWeights` and
  `setPredictor` move no forest's sd. A control taken from a sampler and given to a fit of a response twenty
  times as wide carries no anchor: the new sampler's is its own (67.71 against 3.386). A warm start refuses
  several forests, and a model of several forests refuses test predictors, so no test row can reach a basis.
- The prior-only path. `sampleTreesFromPrior` with `sampleLeafParametersFromPrior` draws each forest's own
  fit at the sd the reader reports (three forests over 3000 draws: 0.997, 1.027, 0.987 of `k.scale`) and
  leaves every coefficient where it was. [`samplePriorPredictive`](../../R/dbarts.R) on a sampler of several
  forests stops inside `predict`, on push 3's build too, in a text about test columns.
- A state. Its glue block is K, the widths, every coefficient, held ones included, and one variance per
  forest ([`readAmplitudeGlue`](../../src/R_interface_bartcore.cpp)). An install takes a forest's
  coefficients where the recipient draws them and its variance where the recipient's size is half-Cauchy
  ([`restoreGlue`](../../src/bartcore/combiner.hpp)); a held block in the state is passed by, the rule
  ["fixed amplitudes"](../../inst/tinytest/test-state-not-model.R) pins. So a state from a drawn sampler goes
  into one that holds (0, 1), the trees arriving and the coefficients staying (0, 1), "exact" `TRUE`. A state
  holds no sd, no scale, no row norm and no kind. (It records its family at the tip, a record dec-B283
  removes; nothing here reads it.)
- The engine, as forest-sd-unit leaves it. [`ForestSpec`](../../src/bartcore/combiner.hpp) states an sd
  or none, a hold and, where `amplitude.prior.variance` is stated, a variance; one function of the engine
  resolves the tip's law from that and from whether the forest has a basis, and the forests' record
  holds four numbers a forest. The K-forest constructor of [`Chain`](../../src/bartcore/chain.hpp) still
  divides by [`basisRowNorm`](../../src/bartcore/chain.hpp)
  ([`mapLeafScale`](../../src/bartcore/chain.hpp)). [`ForestAmplitudePrior`](../../src/bartcore/combiner.hpp)
  holds one variance for a forest's whole block, which is all one numeric column or two levels need.
- The shipped header. [`DBARTS_C_API_LIST`](../../inst/include/dbarts/dbarts.h) creates no sampler from a
  forest's statement and reads no forest's spread, scale or hold; a per-draw callback sees the coefficients
  as `glue`. Nothing in this slice changes a signature there.
- The reader, as forest-sd-unit leaves it. `$getLeafPrior(f)` on a forest of several gives `leaf.prior`
  (the statement), `sd`, `sd.stated`, `prior.sd.of` ("amplitude scale" or "forest total"), `k.scale`,
  and `amplitude.prior.scale` or `amplitude.prior.variance`, `leaf.scale.factor`, `leaf.scale.divisor`,
  `basis.row.norm`; with no forest, a list named by label. `extract(type = "leaf.prior.sd")` is the list of
  the `sd` entries by label, and the entry itself for one forest. `print` of a fit, `show` of a sampler
  and the verbose summary say nothing of any forest's sd.
- A stand-in for the slice: the forest-prior-args design's prototype of the law, ported to the tip (the
  engine, the bridge and R in the tip's spelling, an sd still in units of L; it also carries the
  kind-by-class rule). On it:
  - The law. 28 rows of the table under "The rule" and of its several-column rows since cut, every kind,
    held and drawn, default and stated, at the second, third and first position, under the three
    families: measured size over the rule's 0.967 to 1.032. tests/cpp builds against it and passes, 350
    lines ok.
  - What keeps its bits. 20 of 20 seeded fits of forests with no basis or a factor, coefficients drawn, are
    `identical()` on the tip and the stand-in (two and three forests, a factor third, stated sds, probit,
    logistic, weights and an offset, `subset`, `bart` with two chains, a logical basis, a swap, a write, a
    copy, a reload, `setResponse`, forest weights and a mask, no forest without a basis); a held forest
    with no basis at its default under probit too. 11 fits with a number or a hold differ. The BCF
    equivalence compare: 14 of 15 scenarios identical, `glue_toggle` not (111 of 209 summaries beyond
    z = 3, the largest 13.95).
  - The suite, at the tip. With its refusals logged and passed by, 15972 results run, one file stops and
    132 fail in 9 files: 72 pin the stored numbers by position (forest-sd-unit's to repair now), 28 the
    reader's entries and labels, 13 the row norm, the coefficient's variance and a start by position, 8
    are dec-A171's refusals, 11 are the kind-by-class slice's. No seeded draw of a drawn forest with no
    basis or a factor moved.
  - Its faults, each a requirement here. A declaration that gives a forest of a data object another basis
    (`dbartsSpec()` over `sampler$data`) keeps the old column's recorded scale: sd 0.317 where the new
    column, a hundred times as wide, has 0.00317. Its engine decides "constant" by a computed standard
    deviation above zero: with the record cleared, `new("dbartsSampler")` on a constant basis refuses 1
    and 0.5 and CREATES 0.1, 61.7 and 1e8 + 0.1, with recorded scales 5e-16, 3e-13 and 5e-7 and own sds
    of 5e15, 7e12 and 5e6 response units; 61.7 + 1e-9 sin(i) is created too. Its one-pass standard
    deviation is far off R's on a column far from zero.
- On a constant column of the suite: push 3's build creates 55 forests on one, 54 of them in
  test-bcf-family.R, and 45 state no sd.
- The exact gate, drafted and run against the prototype's own build. One binary predictor, one tree a forest,
  sigma held; the leaves integrate out and the only quadrature is over the coefficient. Nine arms of one
  numeric column (a 0/1 column and a column of standard deviation 0.107, default and stated, drawn and
  held; a column with a fifth of its rows at weight zero, drawn and held). With the verdict the previous
  planner then tried, fail when the gap exceeds 0.0005 plus 4 standard errors, at 32 seeds of 25000 kept
  sweeps: the arms pass (worst 0.52 of the bound); a recorded scale times 1.03 fails every arm tried
  (4.6 to 10.6), times 1.05 and 1.10 fail everywhere; n for n - 1 in the standard deviation, a scale
  times 0.987, fails the drawn arms only just (4.07, 4.06) and passes the held ones (2.6, 1.5). Under
  the first plan's fixed tolerance of 0.01 a scale 3 percent off passed every arm and 5 percent off
  passed the held ones (the critic's run). Against an oracle written to a wrong law, as gaps: 0.674 left
  on a held coefficient 0.03 at the default, 0.10 and 0.21 at sds of 0.6 and 0.2 of the unit; a
  coefficient variance of 1/2, 0.11; a held coefficient at 0, 0.27; a centred column 0.51; a stated sd
  read per standard deviation 0.14 to 0.65; a default per unit 0.06 to 0.85; the scale taken over the
  rows of positive weight 0.18 and 0.47; a held plain forest that ignores an sd of 0.3 of the unit,
  0.033. An arm of three trees a forest with K = 3 (a forest with no basis held, a numeric column drawn
  or held, a two-level factor held third; the oracle sums over how many of each forest's three trees
  split, 64 combinations, closed form) is drafted and passes (2.35 and -0.02); a leaf's sd not divided
  by sqrt(3) gives 167 and 338, the default taken for K = 2 gives 70 and 52, a scale times 1.03 gives
  9.9 and 12.6; 18 seconds drawn and 8 held. One functional of about a hundred sat at z = -3.07 at 32
  seeds and at -1.49 at 64 seeds of 50000: noise.
- What the BCF gates hold. [bcf-exact.R](../../benchmarks/R/bcf-exact.R),
  [bcf-exact-weak.R](../../benchmarks/R/bcf-exact-weak.R),
  [bcf-exact-restricted.R](../../benchmarks/R/bcf-exact-restricted.R) and
  [bcf-latent-exact.R](../../benchmarks/R/bcf-latent-exact.R) each hold a coefficient in some arm, write
  ["scaleTau"](../../benchmarks/R/bcf-exact.R) with the 0.674 whether it is held or not, and pass
  `amplitude.prior.variance = 0.5`; [sbc.R](../../benchmarks/R/sbc.R) has a
  ["fixedGlue"](../../benchmarks/R/sbc.R) arm. Run on the stand-in, in `quick`, with their oracles as they
  are: the three gaussian gates pass under the new law (largest gaps 0.0034, 0.0224 against 0.03, and
  0.0004), so they cannot see it. Nor can a small sd make them: by bcf-exact.R's own oracle, in its
  design (150 rows, noise 0.3), the two laws differ by at most 0.0042 with sd.control 0.3 and sd.moderate
  0.2, and by 0.0056 at sd.moderate 0.1, against its tolerances of 0.05 and 0.015. The latent gate fails
  in six of the seven arms that hold the treatment forest's coefficients (worst |z| 5.34 to 10.70; the
  seventh, three forests under probit, 3.97 and inside its bound) and passes in the two that hold only
  the prognostic one at its default (1.72, 2.66), in 216 seconds.
- Consumers. bartCause's `bcf()` (dbarts-1.0, 6c1bff9) passes `amplitude.prior.variance = b.prior.variance`
  on its treatment forest and holds either forest with `amplitude = fixed()`; it reads
  `getLeafPrior(1L)$response.scale` and `$response.shift` and the coefficients as `glue`, and nothing else of
  this surface. Two of its test files build the same sampler by hand with `amplitude.prior.variance = 0.5`.
  Its treatment basis is a two-level factor from the kind-by-class slice on. stan4bart (bartcore,
  a9d081b), treatSens (dbarts-1.0, aecec71) and bairrtt (main, 3f57f61) declare one forest with
  `dbartsForests$forest(n.trees = )` and read no coefficient, no forest's sd and no leaf prior of a model
  of several forests (searched again 2026-10-07).

## The rule

s is a forest's `sd`, in the response's units; U is the unit forest-sd-unit records (sd(y - offset) over the
rows kept under gaussian, 1 under probit, pi / sqrt(3) under logistic); d is sqrt(2 / K). Term f at row i is
(B_f(i, .) . a_f) F_f(x_i): the basis row times the forest's coefficients, times the forest's own fit. The
kind of a forest is what forest-kind-by-class records: nothing, the two levels of a factor, or a numeric
column.

| the forest multiplies | coefficients | drawn | held, `amplitude = fixed()` |
|---|---|---|---|
| nothing | one | F has sd U; the size of a is half-Cauchy with median s / U: a F has median size s | a is 1; F has sd s, exactly |
| the two levels of a factor | one a level | each a_l is N(0, 1/2); F has sd s / 0.674: (a_2 - a_1) F has median size s | a is (0, 1); F has sd s, exactly, the second level against the first |
| a numeric column w | one | a is N(0, 1); F has sd s / 0.674: a F has median size s, per unit of w | a is 1; F has sd s, exactly, per unit of w |

1. Stated. `sd = s` is one number, per unit of the column as it was given where the forest multiplies a
   number, and nothing is read from the data.
2. Not stated. No basis: 2 U under gaussian, U under probit and logistic. A factor: d U. A number:
   s = d U / sd(w). Held or drawn alike.
3. The scale, sd(w): the sample standard deviation of the column, n - 1 in the divisor, unweighted, over
   every row the sampler holds when it is first created. Those are the rows left by `subset` and the
   na.action; rows of weight zero count, and no later mask enters. The engine takes it, once, and
   R records it with the data; every later construction is handed the record and takes nothing again. A
   constant column has no scale, and neither has one whose standard deviation is under 2^-26 of its
   largest absolute value: where a default would need one it is refused by name.
4. The column is never centred and never divided: the engine multiplies by it as given, at the fit, after a
   swap and at new rows. A default fit on w and on c w is the same fit, for any c that leaves both
   representable: a row leaves a forest's update when its multiplier is at most 2^-26 of the largest
   absolute multiplier that forest has in the sweep, and for no smaller reason.
5. A held value goes by the forest's kind and never by where the forest stands.
6. The variance of a coefficient is not an argument: `amplitude.prior.variance` goes.
7. A swap moves the basis and nothing else: no scale, no sd in force, no coefficient, no tree.
8. A write, `$setLeafPrior(forests = )`, states s, as forest-sd-unit has it; the record of a scale stays.
9. Where a drawn coefficient starts is not part of the law. A forest with no basis or a factor starts where
   it starts today, which is what keeps its draws; a numeric coefficient starts at 1.

Before and after, per path.

| path | before | after |
|---|---|---|
| creation at every door: a `forest()` term at `bart()` and `dbarts()`; a `forests` list on a formula, on a matrix, in `dbartsSpec()`; `dbartsData(bases = )` | the table of Context | the rule; a numeric forest that states no sd has its scale taken and recorded |
| creation from a data object that carries a scale (`sampler$data`) | | used, for a forest that states no sd; checked there and nowhere else |
| a declaration that gives a forest of such an object another basis | | that forest's record is dropped and taken from the new column |
| `$setForestBasis` | the row norm is taken again: the forest's own sd follows the new column | nothing but the basis moves |
| `$setLeafPrior(forests = )` | the number goes to one of two channels | the forest is stated at s |
| `copy()`, a reload, `new("dbartsSampler", control, model, data)` | built from the column then in force | built from the record: the creation's prior, bit for bit |
| `setState`, a state in a reload | carries no prior | unchanged; a held block in the state is passed by, as today |
| `predict` at new rows, `fitted`, `extract` of draws | coefficients times the basis as given | unchanged |
| `sampleTreesFromPrior`, `sampleLeafParametersFromPrior` | each forest's own sd; coefficients left | unchanged, at the rule's own sd |
| `samplePriorPredictive` | stops inside `predict` | refused by name |
| `setData`; a warm start | refused | refused |
| `setResponse`, `setOffset`, `setWeights`, `setPredictor`, `setActiveRows` | move no forest's sd | unchanged |
| a control or a `forests` list carried to another fit | the forests' record is cleared | unchanged; a scale is on the data and never on the control |
| a state stored before | | installs as any state; it holds nothing of the law |

## The engine and 1.0-0

dec-B281 and dec-B282 leave a 1.0 user three things to multiply a forest by: nothing, the two levels of
a factor, one numeric column. The engine is not R-specific and holds every default; what it keeps
general and what is simply not built:

- Kept as it is: a coefficient for each column of any block, their draw as one block under one variance
  ([`drawForestAmplitude`](../../src/bartcore/combiner.hpp), not edited by this slice), the state's glue
  block for any widths, a swap and a prediction at any width. A two-level factor is two columns there.
- Kept, and new here: the forest's kind is stated to the engine (push 1), and the law is one function
  with a row for each kind. A two-column block is a factor's because it was said to be, never because of
  its width, so a numeric block of two columns, whichever of its meanings is chosen after the release, is
  a new row and no change to how a 1.0 model is read.
- Not narrowed: a drawn levels block of more than two levels keeps the factor's law in the engine, which
  it has today. R refuses it; nothing in tests/cpp is written on it by this slice.
- Not built: a prior scale for each coordinate of a block, a default for each of several columns, a
  stated sd for each column, a held value for several numeric columns or for three or more levels, any
  law for a numeric block of several columns. The engine throws on the last three; R never reaches
  them. dec-B272 has the engine hold one scale a column from this slice on so that no stored state and
  no compiled interface changes when a size for each column is added. Under dec-B282 no 1.0 model has a
  second column to scale, and what that sentence protects holds without the code: a state stores no sd
  and no scale, and the shipped header reads none.

## The reader, `extract` and the printed block

One forest's entry from `$getLeafPrior(f)` after push 2, on a model of several forests. The figures are for a
response with sd 1.452, a dose with sd 0.0204 and an age with mean 50.88 and sd 17.55; in the table each
forest stands beside one forest with no basis (K = 2), and the list and the block below are one model of
four forests, `forest(x1 + x2) + forest(x1, basis = dose) + forest(x1, basis = age) +
forest(x1, basis = factor(z), amplitude = fixed())` (K = 4). Each literal a test pins is computed in the
test from its fixture.

| forest | `leaf.prior` | `sd` | `sd.stated` | `multiplier` | `amplitude` | `prior.sd.of` |
|---|---|---|---|---|---|---|
| no basis, default | `forest()` | `2.904` | `FALSE` | "none" | absent | "forest (prior median)" |
| no basis, held, default | `forest()` | `2.904` | `FALSE` | "none" | `fixed()` | "forest" |
| a numeric column, default | `forest()` | `c(dose = 71.18)` | `FALSE` | "numeric" | absent | "forest, per unit of basis (prior median)" |
| the same, `sd = 30` | `forest(sd = 30)` | `c(dose = 30)` | `TRUE` | "numeric" | absent | the same |
| a numeric column, held, `sd = 30` | `forest(sd = 30)` | `c(dose = 30)` | `TRUE` | "numeric" | `fixed()` | "forest, per unit of basis" |
| factor, default | `forest()` | `1.452` | `FALSE` | "factor" | absent | "level difference (prior median)" |
| factor, held | `forest()` | `1.452` | `FALSE` | "factor" | `fixed()` | "level difference" |

- `leaf.prior`, `sd` and `sd.stated` are forest-sd-unit's. What this slice changes of them: a numeric
  forest's `sd` is one number named as its column is, where push 3 of written-surface names the column,
  and unnamed otherwise; a forest with no basis or a factor keeps one unnamed number. That is
  dec-B275's shape, a vector by column, with one column.
- `basis.scale`: on a numeric forest whose data carries one, the recorded standard deviation, named as
  `sd` is.
- New: `multiplier` and `amplitude` (push 1). Kept as they are: `leaf.model`, `prior.mean`, `k.scale` (the
  forest's own sd), `response.scale`, `response.shift`. Gone: `amplitude.prior.variance` (push 1);
  `amplitude.prior.scale`, `leaf.scale.factor`, `leaf.scale.divisor`, `basis.row.norm` (push 2).
- `extract(fit, type = "leaf.prior.sd")` on that fit of four forests:
  `list(forest1 = 2.904, dose = c(dose = 50.33), age = c(age = 0.0585), "factor(z)" = 1.027)`; with
  `forest = "dose"` it is `c(dose = 50.33)`, the entry itself.
- The printed block, in `print` of a fit, `show` of a sampler and the verbose summary (push 3):

      forests, sd in the response's units (its sd 1.452):
        1 forest1: no basis; sd 2.904 (default)
        2 dose: numeric basis, mean 0.01994, sd 0.0204, range 0 to 0.0698; sd 50.33 per unit, 1.027 per sd (default)
        3 age: numeric basis, mean 50.88, sd 17.55, range 19 to 88.4; sd 0.0585 per unit, 1.027 per sd (default)
        4 factor(z): factor basis, 2 levels; sd 1.027 between the two levels (default); coefficient held

  A stated sd reads `sd 30 per unit, 0.612 per sd (stated)`. Under probit the header is `forests, sd on the
  latent scale (error sd 1):` and under logistic `(error sd 1.814)`. A column's mean, sd and range are
  those of the column in force; where a default was taken from another column the line ends `(default,
  taken from a column of sd 0.0204)`. A mean is printed through `zapsmall` against the sd. This is the
  one sign a 0/1 number (dec-B263: at 20 percent treated beside one other forest, `mean 0.1967, sd
  0.3981, range 0 to 1; sd 3.647 per unit` where the factor's line says 1.452), a column far from zero
  (dec-B264: `mean 50.88, sd 17.55`), one wild value (the range) and a stated sd in the wrong units get:
  nothing warns. For a default, "per sd" is the same number on every numeric line of a model; it is
  there for the stated case.

## Refused forms, with their texts

Base R's style, as in written-surface. `<f>` is `forest 2`, or `forest 2 ("dose")` where the forest has a
label. The push that adds each is in brackets.

    [1] 'amplitude = fixed(2)': a held coefficient is 1 for a forest with no basis, and 0 and 1 for the two levels of a factor; fixed() takes no other value. Write fixed(), and state the forest's size with 'sd'
    [2] 'amplitude = fixed(2)': a held coefficient is 1 for a forest with no basis and for a numeric column, and 0 and 1 for the two levels of a factor; fixed() takes no other value. Write fixed(), and state the forest's size with 'sd'
    [2] <f> states no sd and its basis is constant (61.7 on every row the fit keeps), so there is no standard deviation to take a default from; state one with sd =
    [2] <f> states no sd and its basis is constant to rounding (its standard deviation is 1e-11 of its largest value), so there is no standard deviation to take a default from; state one with sd =
    [2] samplePriorPredictive does not support a sampler of several forests: their coefficients are not drawn from the prior here
    [2] 'basis.scale' entry 2 must hold one positive finite number for each column of that forest's basis: it has 2 and the basis 1 column

An argument that does not exist gets R's own error: `unused argument (amplitude.prior.variance = 0.5)` from
push 1, and `updateBasisScale =` on `$setForestBasis` as today. Behind these, from the engine and the bridge:
"a held coefficient block must be a forest's with no basis, a two-level block or one numeric column"; "a
forest with no basis takes none"; "a basis column is constant, so no default sd can be taken from its
standard deviation"; "a numeric basis of several columns has no law"; "the sampler refuses this basis for
the forest", where today the bridge drops the engine's answer. Gone with push 1:
forest-defaults-by-kind's two interim texts, for a held forest with no basis that is the second and for
a held basis of two columns that is not;
"'amplitude.prior.variance' is the prior on a basis forest's amplitudes". Gone with push 2: dec-A171's
refusal of one numeric column. Kept as they are: every text of the kind-by-class slice, whose refusals
of a basis's size come before any of these.

## The help's `sd` item, as push 2 leaves it

    sd: The size of what this forest contributes, in the response's units: the units of y for a
    continuous response, and of the latent index under "probit" and "logistic", where the link's own
    error has standard deviation 1 and pi / sqrt(3). One number; a name on it is ignored. What it is
    the size of depends on what multiplies the forest.

    No basis. The forest enters as a f(x), and sd is the size of a f(x).
    A factor (a factor, a character vector or a logical vector, of two levels). The forest enters with
    one coefficient for each level, and sd is the size of the difference between the two: the treatment
    effect where the levels are control and treated.
    A number w. A column that is an ordinary covariate should be centred: basis = scale(w), or
    basis = I(w - c) for a reference value c. The other forests describe the response where w is 0, and
    a column that never comes near 0, an age or a year, makes them an extrapolation and the fit poor,
    with nothing in the chain to show it. A dose or an exposure, where 0 means none, is left as it is.
    The forest enters as w (a f(x)), and sd is the size of a f(x), the change in the response for one
    unit of w; the term at an observation has size |w| sd. The column is used as given: it is never
    centred and never divided, so a stated sd is per unit of w, and the same sd with w in grams and in
    milligrams states two priors a thousand times apart.

    Where the coefficient is held (amplitude = fixed()) "size" is the prior standard deviation, exactly.
    Where it is drawn the size is itself uncertain, and sd is its prior median.

    Not stated, sd is 2 sd(y) for a forest with no basis under a continuous response and one unit of the
    link's error otherwise, whether its coefficient is drawn or held; sqrt(2 / K) of sd(y), or of that
    unit, for a factor in a model of K forests; and for a number that same size for one standard
    deviation of the column, sqrt(2 / K) sd(y) / sd(w) per unit. A default fit is then the same fit
    whatever units w is in. sd(w) is taken once, from the rows the sampler is created with, rows of
    weight zero included, and kept with the data: replacing the column later does not move the prior. A
    caller who replaces a column and wants the default of the new one states it,
        p <- sampler$getLeafPrior(2)
        sampler$setLeafPrior(forests = list(forest(), forest(sd = p$sd * p$basis.scale / sd(w.new))))
    A constant column has no standard deviation and needs sd stated.

    A treatment belongs in a factor or a logical vector, where sd is the treatment effect. Written as a
    0/1 number it is scaled as any number: sd(z) is 0.5, 0.4 and 0.22 at 50, 20 and 5 percent treated, so
    the default effect is 2, 2.5 and 4.6 times the one factor(z) gets;
    sd = sqrt(2 / K) * sd(y) on the number states the factor's.

    A held coefficient cannot adapt to an sd that is too small. At the default it is sized for a slope of
    about one sd(y) per sd(w); where the effect may be larger, under "probit" and "logistic" above all,
    state sd or let the coefficient be drawn.

After push 1 the item has the paragraphs on a held coefficient and on a forest with no basis and a
factor, and says of a number what forest-sd-unit's help says. Until push 3 it does not say that `print`
shows what each forest resolved to; push 3 adds the sentence.

## Constraints

- A fit whose forests have no basis or a factor, coefficients drawn, draws what it draws on the base build,
  to the bit, in every push, stated sds included. In push 1 a fit with a drawn numeric multiplier does too.
- The help promises a default fit that does not depend on a column's units only on a tip that has the
  relative bound: the bound is push 2's first commit and the sentence its last.
- The law is resolved in one function of the engine and nowhere in R: R states what the caller stated and
  the kind, and reports what the engine returns. The engine holds every default, the column's scale
  among them: R computes no standard deviation that the engine uses.
- A forest's law is resolved once for each construction, by the chain, and the combiner is handed the
  result; no second derivation can differ from the first.
- "Constant" is decided by comparing values, in R and in the engine, never by a computed standard
  deviation being zero.
- [`inst/include/dbarts/dbarts.h`](../../inst/include/dbarts/dbarts.h) is not edited and
  [`DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h) does not move.
- The stored state is not edited: no block, no name, no encoding;
  [`stateFormatVersion`](../../src/R_interface_bartcore.cpp) stays.
  [`serializeGlue`](../../src/bartcore/combiner.hpp), [`restoreGlue`](../../src/bartcore/combiner.hpp) and
  [`glueIsValid`](../../src/bartcore/combiner.hpp) are not touched.
- `data@bases[[f]]` is not touched: no column is centred, divided or given an attribute.
- The unit, the conversion of a stated sd, the statement the engine is handed, the reader's `leaf.prior`,
  `sd` and `sd.stated` and `extract`'s shape are forest-sd-unit's and are carried, not rebuilt. A default
  never passes through the unit on its way in.
- The kind is read from the one function forest-kind-by-class adds, and its swap table stands as it is:
  no forest changes width, so this slice adds no row to it.
- [`drawForestAmplitude`](../../src/bartcore/combiner.hpp) is not edited, and the lines
  ["m13"](../../benchmarks/R/mutation-battery.R) and m14 quote in it stay as written; the entries that
  quote lines this slice rewrites move with them.
- No refusal is written, and no reading kept, for a sampler or a fit that only an earlier build of this
  branch could have saved; the bridge's check of the record's length stands.
- Nothing of `updateBasisScale`, of `leaf.prior` on `forest()`, of a heavy-tailed coefficient or of
  several columns rides along.
- Each push leaves the help saying what the code does. Base R calls stay within DESCRIPTION's R floor.

## Pushes

Three, each gated on its own and each a coherent tip.

1. Held coefficients. The engine is told each forest's kind; the held values go by kind, a held forest's
   sd is exact, forest-defaults-by-kind's two interim refusals are lifted, and
   `amplitude.prior.variance` goes. The numeric law is untouched, so dec-A171's refusal of a held
   numeric column stays. The exact gate's script is created with its held arms. bartCause's edit the
   same day.
2. The numeric multiplier. First commit, gated alone and bit for bit: the bound under which a row leaves
   a forest's update becomes relative. Then the law: the column as given, the normal coefficient, the
   default for one standard deviation of the column with its record, a held numeric column, the
   reader's and `extract`'s last change. dec-A171's refusal is lifted. The gate's numeric arms and four
   baseline scenarios. The design note is completed here.
3. The printed block. R only.

No tip between them fits a model nobody asked for: after push 1 a held numeric column is still refused,
not held under half a law; push 2 lands the numeric law, its reader and its help together. Push 2's tip
lacks the printed block for as long as push 3 takes; its help does not promise one.

## Steps

"Fails today" is what the base build does where the test expects otherwise. New names are for the reader's
sake. Calls are in push 3 of written-surface's spelling; sds are in the response's units. Fixture, unless
said: 150 rows with x1, x2, a 0/1 z, its factor zf, dose (sd 0.02), age (mean 50, sd 17), an offset column
and weights of which ten are 0; a response whose standard deviation is neither 1 nor its range.

### Push 1: held coefficients

1.1 The engine. [`ForestSpec`](../../src/bartcore/combiner.hpp) states a kind (nothing, levels,
    numbers) beside what forest-sd-unit has it state, and
    [`expandForestSpecs`](../../src/bartcore/combiner.hpp) states BCF's two forests as nothing and
    levels. The law's function (`forestLaw`) takes the held rows of "The rule" for a forest with no basis
    and a two-level block: a held forest with no basis has own sd s, by the leaf scale and not the
    half-Cauchy median; a held two-level block has no divisor. A held block is set to its kind's value
    when the combiner is built, at any position, after
    [`rebuildAmplitudeLayout`](../../src/bartcore/combiner.hpp), which keeps the drawn starts; a held
    block with no value throws there: a levels block of other than two columns, and in this push every
    numeric block. Thrown too: a forest stated to have no basis handed one, and a levels or numeric
    forest handed none. [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp) on a held forest with no
    basis restates the leaf scale. The stated variance leaves
    [`ForestSpec`](../../src/bartcore/combiner.hpp) and
    [`AmplitudeSpec`](../../src/bartcore/combiner.hpp): a level's is 1/2, and a numeric column keeps the
    tip's 1/2 and its row norm until push 2. Every struct here is read by objects that do not track
    headers: `--preclean`. Each fixture of tests/cpp states the kind of every forest it gives a basis;
    a comparison that holds no coefficient is unchanged by that, to the bit.
    Tests, tests/cpp (`testHeldCoefficients`): three forests (none, a two-level block, a numeric column
    that is drawn) in each of their six orders, the first two held: the held values by kind at every
    position (fails today in four of the six); each held forest's leaf scale equal to its sd over the
    unit times the anchor, to 4 ulp, with sds that differ from each other, from the unit and from 1
    (fails today: the divisor, and a forest with no basis at the anchor); the drawn numeric forest's
    leaf scale the base build's, to the bit; the throws for a held three-level block and a held numeric
    block, and the two for a kind that does not fit its basis; a swap of a held two-level block to
    another of its width leaves (0, 1); a state stored from a twin that draws installs and leaves
    (0, 1); a write to a held forest with no basis gives the leaf scale of one created at that sd.
1.2 The bridge and the kind. [`parseData`](../../src/R_interface_bartcore.cpp) gives each forest its kind
    from the data object: no basis; a basis with the levels forest-kind-by-class records; a basis
    without. [`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp) no longer reads a variance.
    [`bartcore_getLeafPrior`](../../src/R_interface_bartcore.cpp) reports the kind the engine holds and
    the hold. Tests, a new file test-forest-held.R, first block: each of six objects (a factor, an
    ordered factor, a character and a logical vector, a 0/1 number, dose) at each of the seven creation
    doors of forest-kind-by-class, and again after `copy()`, a reload and
    `new("dbartsSampler", control, model, data)`: the kind the bridge reports is the kind that slice's
    function gives for the data object (a kind dropped between the data and the engine reads "numeric"
    for a factor at that door).
1.3 R. [`forest`](../../R/model.R) loses `amplitude.prior.variance`, and
    [`validateForestKnobs`](../../R/model.R), [`resolveForests`](../../R/model.R) and the multi-forest
    block of [`resolveSamplerSpec`](../../R/spec.R) their lines for it. forest-defaults-by-kind's
    function for the held shapes keeps one row, dec-A171's, and loses the two by position: a forest with
    no basis and a two-level factor are held wherever they stand.
    [`validateForestAmplitude`](../../R/model.R) takes the [1] text. The reader
    ([`reportLeafPrior`](../../R/dbarts.R)) gains `multiplier` and `amplitude` (`fixed()` or absent),
    takes the four `prior.sd.of` values of a forest with no basis or a factor, and loses
    `amplitude.prior.variance`.
    Tests, test-forest-held.R:
    - By kind, at every door and position. `amplitude = fixed()` on a forest with no basis, on
      `factor(z)` and on a logical, in a formula term, a `forests` list on a formula and on a matrix,
      `dbartsSpec()` and a data object's `bases`; as the first, second and third forest, and on a forest
      with no basis written second: `getForestAmplitudes()` is 1, or (0, 1), before and after 50 sweeps
      (fails today: refused at the positions forest-defaults-by-kind closed; before that slice (1, 1)
      third and first, and 0 for a plain forest second).
    - Exact. The reader's `k.scale` of a held forest equals its `sd`: stated 0.3 and default on a forest
      with no basis, stated 0.6 and default on a two-level factor, under each family, to 1e-12 (fails
      today: the response's scale, and s / 0.674). The default of a held forest with no basis is 2 U
      under gaussian, the drawn default's number.
    - The size, from the engine's prior draws. Always on, one a family, 1500 sweeps and 10 percent: the
      second level's term of a held factor and a held forest with no basis have the stated sd. `at_home`,
      5000 sweeps a row and 4 percent: the same at three sds.
    - Refused: `amplitude.prior.variance` is R's unused-argument error at both doors and through
      `do.call`; `fixed(2)` with its text.
    - Not moved: a model with a numeric column drawn beside a held factor keeps the base build's row norm
      (`k.scale` of the numeric forest to the bit), and dec-A171's refusal stands at creation.
1.4 The gates.
    - A new script, benchmarks/R/forest-law-exact.R, which push 2 extends: one binary predictor, sigma
      held, the leaves integrated out in closed form; no quadrature in this push's arms. Arms: a forest
      with no basis held at sds of 0.3 and 2 of the unit; a two-level factor held at 0.6, at 0.2 and at
      its default; both held; and one arm of three trees a forest with K = 3, a forest with no basis
      held, a two-level factor held second and another held third, the oracle summing over how many of
      each forest's trees split. Each compares the posterior mean of each forest's term in each cell and
      the probability that each tree splits with the closed form written from "The rule"; the oracle
      reads nothing off the sampler but the unit and the shift. The verdict goes by z: an arm fails
      when a gap exceeds 0.0005 plus 4 standard errors over seeds; `quick` is 32 seeds of 25000 kept
      sweeps, full 64 of 50000. With about a hundred functionals a run, a lone z between 4 and 5 is
      expected in a few runs of a hundred: the script reruns that arm once on a second set of seeds and
      fails if the same functional is out again. Added to
      [exact-gates.yaml](../../.github/workflows/exact-gates.yaml)'s list.
    - The four BCF gates: the `amplitude.prior.variance` argument goes, and each oracle takes the
      treatment forest's sd with no 0.674 in an arm that holds its coefficients and the prognostic
      forest's sd as stated in an arm that holds a, so that each says what the engine does. The three
      gaussian ones are not this push's oracle: they pass under either law (Context) and cannot be made
      to see it by a small sd. [bcf-latent-exact.R](../../benchmarks/R/bcf-latent-exact.R) is, under
      probit and logistic: run it once on the push's build with its oracle left as it is, six arms
      failing, and then restated, none.
    - [sbc.R](../../benchmarks/R/sbc.R): the argument goes, and its fixed-glue arm draws its truth at the
      held scales; the arm constructs and runs its smallest setting.
1.5 The baseline. Run the BCF compare against `bcf-equivalence-1b7d730c.rds`: 14 scenarios identical, with
    no `max |z|` line, and ["glue_toggle"](../../benchmarks/R/bcf-equivalence.R) not. Re-record; the
    MANIFEST row names the oracle (rule P17), which is an identity and a gate: the re-recorded draws equal,
    to 1e-10, those the base build gives the same scenario with its treatment forest stating 0.674 of the
    default (pair row 23's identity, run for the scenario itself); and the held-factor arms of
    forest-law-exact.R and of bcf-latent-exact.R pass, with their gaps.
1.6 Respell. `amplitude.prior.variance` on the test lines that state it (31 lines of 6 files at the tip,
    recounted after the kind-by-class slice's removals): where it states 0.5 the argument is dropped and
    the fit is the same; where it states another value the test pinned the tip's law and goes with it,
    one pin of R's error taking its place. Five benchmark scripts drop the argument.
    forest-defaults-by-kind's pins of its two interim texts by position become step 1.3's creations. The
    reader's pins of `prior.sd.of` and `amplitude.prior.variance` for a forest with no basis or a factor
    take the new entries. Run the suite first and repair what it shows.
1.7 Help and records. man/forest.Rd: the usage, the `amplitude.prior.variance` item removed, the
    `amplitude` item (the held values, two in this push, wherever the forest stands; that a held
    coefficient makes `sd` exact; that a held forest with no basis has the drawn default, 2 sd(y) under a
    continuous response), the `sd` item's sentences for a held forest, the Details paragraph on the
    budget. [`dbartsSampler$getLeafPrior`](../../man/dbartsSampler-class.Rd) with its docstring. A new
    docs/design/forest-multiplier-law.md with its index row: the rule's held half, the changed sequences
    with their oracles, the re-recorded scenario, what the engine keeps general. docs/design/bcf.md and
    multiplier-combiner.md where they give the 0.674 to a held forest or name the variance. TODO:
    `forest-prior-args`.
1.8 bartCause, same day (its own commit on dbarts-1.0; its basis has been `factor(as.integer(z),
    levels = 0:1)` since the kind-by-class slice). In R/bcf.R: `b.prior.variance = 0.5` leaves the
    formals of `fitBCF` and of `bcf`; `amplitude.prior.variance = b.prior.variance` leaves the treatment
    forest's call; `b.prior.variance = b.prior.variance` leaves the call of `fitBCF`; and beside its
    refusal of `forests` and `bases`, `bcf` refuses the argument by name ("'b.prior.variance' is no
    longer an argument: the treatment coefficients' prior variance is fixed; state the effect's size with
    'sd.moderate'"), since its dots would otherwise take it in silence. man/bcf.Rd: the usage line and
    the item go; `update.a, update.b` reads "when FALSE the coefficient is held, at 1 for the prognostic
    forest and at 0 and 1 for control and treated, and sd.control or sd.moderate is then that forest's
    prior standard deviation, exactly". tests/testthat/test-14-bcf.R and test-03-responseFit.R: the two
    hand-built samplers drop `amplitude.prior.variance = 0.5`, the same fit. A `bcf()` fit with
    `update.a = FALSE` or `update.b = FALSE` moves; no test there pins a draw under either. stan4bart,
    treatSens and bairrtt are not edited, at this push or the next two.
1.9 Mutations (Verification).

### Push 2: the numeric multiplier

2.1 First commit: the bound. In
    [`AmplitudeForestCombiner::formForestResponse`](../../src/bartcore/combiner.hpp) one pass takes the
    largest absolute multiplier the forest has over the rows, and a row leaves with zero weight and zero
    response when its multiplier is at most
    [`zeroMultiplierTolerance`](../../src/bartcore/combiner.hpp) times that; a forest whose every
    multiplier is zero loses every row, as today. The comment above it says what the bound now caps:
    the division amplifies by at most 2^26 relative to the forest's largest multiplier. No other
    combiner's method is edited. Tests, tests/cpp (`testRelativeRowBound`): a numeric column w against
    2^-40 w, drawn, 50 sweeps: identical draws (fails today: blind from 2^-30); a held (0, 1) block
    still gives its first level's rows zero weight; a forest whose coefficient is exactly zero gives
    every row zero weight and finite sums; a row whose multiplier is 2^-30 of the largest leaves, and
    one at 2^-20 stays. tinytest, always on: dose against 2^-40 dose, drawn, under the tip's law, whose
    row norm scales with the column: identical draws. (Run so far on the stand-in's law only. If the
    tip's row norm does not scale to the bit, this commit's twins assert which rows keep their weight,
    and the bitwise twins are step 2.3's and 2.8's.) Gated alone before anything else of the push is
    written on top: tests/cpp; the suite; on a reference build the four snapshot files and the three
    compares, every scenario identical, counted; `bench-sampler.R compare` against the commit before it,
    the one hot path this slice touches (one more pass over the rows each time a forest's response is
    formed).
2.2 The scale. One function of the engine (`basisColumnScale`) gives a column's sample standard deviation
    over every row, n - 1, the mean taken first and the deviations from it squared, as R's `sd()` does.
    It throws for a column with no scale: fewer than two rows; every value equal to the first; a
    standard deviation at most 2^-26 of the largest absolute value. Tests, tests/cpp: columns of sds
    0.25 and 8 give exactly those; a 0/1 column at 20 percent; a column with half its rows zero equals
    the literal over ALL rows; 1e7 + i for i in 1 to 150 equals the literal R gives to 1e-10 relative
    (fails with one pass); the constants 1, 0.5, 0.1, 61.7 and 1e8 + 0.1 each throw (the last three fail
    with a test of the computed sd against zero); 61.7 + 1e-9 sin(i) throws; one row throws.
2.3 The law. [`ForestSpec`](../../src/bartcore/combiner.hpp) gains a scale (`basisScale`, not a number
    for "take it"). `forestLaw`'s numeric row is "The rule"'s: s stated, or the default over the scale
    handed or taken; the forest's own sd s / 0.674 where the coefficient is drawn, under variance 1 and
    starting at 1, and s where it is held at 1. A numeric block of several columns throws, and a
    fixture of tests/cpp that builds one is restated on one column, or as a levels block where that is
    what it stood for.
    [`basisRowNorm`](../../src/bartcore/chain.hpp) goes with everything that read it. One function of
    the chain (`applyForestLaw`) sets a forest's leaf scale from the sd in force, and the constructor
    and [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp) both call it.
    [`Chain::setForestBasis`](../../src/bartcore/chain.hpp) installs the basis and derives nothing.
    [`ForestCalibration`](../../src/bartcore/chain.hpp) reports the sd in force, the scale held and the
    stated flag, and no longer a row norm, a factor or a divisor.
    Tests, tests/cpp (`testNumericMultiplierLaw`):
    - Literals. K = 3 under each family, drawn and held: a forest with no basis, a two-level block and a
      numeric column of sd 0.25, each forest's own sd, its divisor, its coefficients' prior variance and
      start, to 4 ulp (the powers of two make the defaults exact).
    - The twin. w against 1024 w and against 2^-40 w, default: identical draws, drawn and held, gaussian
      and probit.
    - Default against stated. A default forest, its scale taken, and a twin stated at the number that
      default resolves to: identical draws over 50 sweeps.
    - The coefficient's prior. With the forest's structure frozen, the coefficient's conditional mean
      and variance over 4000 draws against the closed form whose prior precision is 1 (fails with the
      tip's 2).
    - Paths in the engine. A chain made again with the first one's scale handed back has its leaf scale
      to the bit, on a basis since swapped to ten times the column; made again with none handed, it
      differs. After a swap, and after a write, the leaf scale times its divisor equals the sd the
      reader reports, to 4 ulp. The throws.
2.4 The bridge. [`parseData`](../../src/R_interface_bartcore.cpp) reads `data@basis.scale` through the
    guard an absent slot needs; [`applyForestBases`](../../src/R_interface_bartcore.cpp) hands a numeric
    forest that states no sd its entry, checked there and only there (one number, finite, positive).
    [`bartcore_getLeafPrior`](../../src/R_interface_bartcore.cpp) carries the scale held;
    [`bartcore_setForestBasis`](../../src/R_interface_bartcore.cpp) raises the engine's refusal.
2.5 The record. A slot `basis.scale` on [`dbartsData`](../../R/A_class.R), beside the levels
    forest-kind-by-class adds: `NULL`, or a list with one entry a forest, `NULL` or a numeric vector with
    one number for each column of that forest's basis, which is one number; read through an accessor
    that answers `NULL` for an object saved without the slot ([`dataRowNames`](../../R/data.R)'s
    pattern). Written in the sampler's `initialize`, where the anchor is recorded, from what the engine
    reports, for a forest whose entry the object lacks; dropped for a forest whose basis a declaration
    replaces on a data object ([`resolveSamplerSpec`](../../R/spec.R)). A write of an sd leaves it.
    Tests, a new file test-forest-scale-record.R. Each path asserts the recorded scale, the reader's `sd`
    and the forest's `k.scale * 0.674` against literals written in the test, with a dose whose sd is not
    1 and two chains:
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
    - Carried and wrong. A record of two numbers, or with a zero, on a forest that states no sd: refused
      with its text; on a forest that states one: created.
    - Stated, then put back. A default forest written `sd = 0.5`, copied and reloaded: `sd` is 0.5 and
      `sd.stated` `TRUE` at each, and the record's entry is still the creation's.
    - Not a path. `predict` at new rows equals the coefficient times the new basis times the forests'
      fits, by `type = "forest"`, to 1e-12, with a new dose a thousand times as wide; a state from a
      sampler created on another column installs and the recipient's record and `sd` are its own;
      `setResponse`, `setOffset`, `setWeights`, `setPredictor` and `setActiveRows` leave the record and
      `sd` `identical()`.
2.6 Refusals and what is accepted. Where the bases are in hand after `subset`
    ([`resolveSamplerSpec`](../../R/spec.R)): a numeric column that is constant, by comparing its
    values, or whose standard deviation is at most 2^-26 of its largest absolute value, on a forest that
    states no sd, by the two [2] texts. [`samplePriorPredictive`](../../R/dbarts.R) refuses several
    forests at its door. dec-A171's refusal goes, and
    [`refuseHeldOneColumn`](../../R/model.R) with it: a numeric column is held at 1.
    Tests: `rep(1, n)`, 0.1, 61.7 and 1e8 + 0.1 on every row, and a column constant on the rows `subset`
    keeps, each refused at four doors and created with `sd = 2`; the same five through
    `new("dbartsSampler", control, model, data)` with the record cleared, where R's check does not
    stand and the engine's must (fails on a check of the computed sd: the last three are created);
    61.7 + 1e-9 sin(i) refused, and 1e7 + i created with the scale R's `sd()` gives; the tests of
    dec-A171's refusal turned over: `amplitude = fixed()` on dose is created at every door, held at 1 as
    the second, third and first forest, with `k.scale` equal to the sd (stated, and
    `sqrt(2 / K) * U / sd(dose)`); swapped to another column it is still held at 1 with the same
    `k.scale`; a held 0/1 number and the held factor of the same column, at one stated sd, have
    `identical()` train draws. `$setForestBasis(2, dose, updateBasisScale = TRUE)` is R's
    unused-argument error.
2.7 The reader, the writer and `extract`. [`reportLeafPrior`](../../R/dbarts.R) names a numeric forest's
    `sd` by its column, adds `basis.scale`, takes the two `prior.sd.of` values of a numeric forest and
    loses the four entries of the tip's channels; [`extractParameter`](../../R/generics.R) is not
    edited and returns what the reader stored; [`leafPriorIsDrawn`](../../R/bart.R) and
    [`defaultLeafScaleVars`](../../R/diagnostics.R) read `multiplier` where they read `basis.row.norm`.
    Tests, a new file test-forest-leaf-prior-reader.R, pinning `sd`, `sd.stated`, `multiplier`,
    `amplitude`, `prior.sd.of` and the round trip, never the class of `leaf.prior`:
    - The literals of the table, for a `bart` fit of the four forests and for the sampler, `sd.stated`
      before and after a write.
    - Through the engine. For every forest of that model and of its held twin, after creation, a write,
      a swap, a copy and a reload: `k.scale`, times 0.674 where the coefficient is drawn, equals `sd` to
      1e-12; for a drawn forest with no basis the half-Cauchy median times the unit does.
    - Read then write, as forest-sd-unit pins it, with a numeric forest among the forests: a default
      numeric forest written `forest(sd = p$sd)`, the name on the number left as it is, has `sd.stated`
      `TRUE` and the next 20 sweeps `identical()` to an untouched twin's.
    - `extract`: `identical()` to the list of the kept sampler's `sd` entries; by a label and by a
      position the entry itself, a numeric forest's with its column's name; a fit with no kept sampler
      the same from what it stored.
2.8 The gates and the tests that run always.
    - forest-law-exact.R gains the numeric arms, nine of one column (Context) and the arm of three trees
      with a numeric column, drawn and held, in the place of its second forest; the oracle written from
      "The rule" with R's `sd()` over every row the sampler holds. Functionals: the posterior mean of
      the term in each cell at each value of the basis; the probability that the multiplied tree
      splits; the probability that the coefficient is within 0.5, 1 and 2. The verdict and the seeds are
      push 1's. What it does not see is said in its header: n for n - 1 in the scale on a held arm
      (step 2.2's literals and step 2.5's 1e-12 hold that), a numeric multiplier under logistic, and
      `subset`.
    - The prior from the engine, tinytest. Always on, one a family, 1500 sweeps and 10 percent: a
      default numeric forest has size `sqrt(2 / K) * U` per standard deviation of its column. `at_home`,
      5000 sweeps a row and 6 percent: default and stated, a column times 1000, an age, a 0/1 number at
      20 percent, held, under each family.
    - The units twin, tinytest, always on: dose against 1024 dose and against 2^-40 dose, default, drawn
      and held, gaussian and probit: identical draws.
2.9 The baseline. Four scenarios are added to [bcf-equivalence.R](../../benchmarks/R/bcf-equivalence.R)
    and recorded: a numeric column at its default; a numeric column through a swap, a write and a
    copy; a held numeric column; a held forest with no basis under gaussian. The 15 in force are
    bitwise on the push's build first, counted; the file then holds 19. The MANIFEST row names the
    oracle of the four: forest-law-exact.R's arms with their gaps, and the units twin. They are there
    for what comes after: the slices that move the forests' record and fold `basis.scale` into the
    data's other records are gated bit for bit, and until now no recorded draw would move if one of
    them took a scale twice or dropped a record.
2.10 Respell and repair. Run the suite first. Expected: the assertions that pin the row norm and a
    numeric coefficient's variance (["medianRowNorm"](../../inst/tinytest/test-bcf-family.R), 13 at the
    tip, fewer after the kind-by-class slice's removals) are rewritten as pins of "The rule" or go with
    it; the 45 creations of test-bcf-family.R on a constant column that state no sd state one; the pins
    of the reader's old entries and labels
    (["forest total"](../../inst/tinytest/test-calibration-midchain.R)) take the new ones. Each test
    that only needs some numeric forest is left, and now fits the law.
2.11 Help and records. man/forest.Rd: the `sd` item above, `amplitude` with a numeric column, `basis`
    where it speaks of a constant column. man/dbartsData.Rd: the slot.
    [`dbartsSampler$setForestBasis`](../../man/dbartsSampler-class.Rd),
    [`dbartsSampler$getLeafPrior`](../../man/dbartsSampler-class.Rd) and
    [`dbartsSampler$setLeafPrior`](../../man/dbartsSampler-class.Rd) with their docstrings: a swap moves
    no prior; the entries. The `extract` page's `leaf.prior.sd` item: a numeric forest's number carries
    its column's name. docs/design/forest-multiplier-law.md completed: the rule, the scale and its
    record on every path, the relative bound and what it caps, what a state and a saved sampler hold,
    the changed sequences with their oracles, and the additions left with the reason each is one.
    docs/design/multiplier-combiner.md rewritten around the rule and the bound; nameable-calibration.md,
    public-surface.md, state-not-model.md and docs/architecture.md where they describe the row norm or
    the data object's slots; the cites of symbols that go, `basisRowNorm` among them, marked `retired:`.
    TODO: `forest-prior-args`, and the entries of "Out of scope".
2.12 Mutations (Verification).

### Push 3: the printed block

3.1 One function (`forestPriorLines`) builds the block from the reader's entries and, for a numeric
    forest, the mean, standard deviation and range of its column as the fit holds it: a sampler's
    `data@bases`, a fit's stored `bases` ([`packageBartResults`](../../R/bart.R)). Nothing more is
    stored on a fit. [`fitSynopsis`](../../R/generics.R) prints it after the tree counts
    forest-defaults-by-kind prints, `show` of a sampler prints it, and a sampler created with
    `verbose = TRUE` prints it once. A fit of one forest prints as it does today.
3.2 Tests, a new file test-forest-print.R: the block of "The reader, `extract` and the printed block" by
    exact string for the model of four forests, at `print`, `show` and creation; a stated forest, a held
    one, probit and logistic headers; after a swap to dose / 10 the column's sd and range are the new
    ones and the line ends with the recorded one; a column centred by `scale()` prints mean 0; a column
    with one value a thousand times the rest shows it in its range; every number in the block read back
    from the text equals the reader's, or R's `mean`, `sd` and `range` of the column, to the digits
    printed; nothing is printed for one forest; no warning anywhere, for an age and for a 0/1 number.
3.3 Help: the `sd` item's sentence that `print` shows what each forest resolved to, with the two lines to
    look for; man/bart.Rd where it describes `print`.
3.4 Mutations (Verification).

## Verification

Every push, against the slice's own library (`R CMD INSTALL --preclean -l <lib> .`, `R_LIBS=<lib>` on every
call; check `dbarts:::buildInfo()$mode` and that the install postdates the source), run in series, at most
two cores (`MAKEFLAGS=-j2`, `EQUIVALENCE_CORES=2`):

- `cd tests/cpp && make && ./test_bartcore`, and again built with `OPT="-O2 -g -fsanitize=address,undefined"`
  under `ASAN_OPTIONS=detect_container_overflow=0`; the R-loaded path under the address sanitizer for the
  push's new test files ([Gate hygiene](README.md#gate-hygiene) gives the commands). Push 3 touches nothing
  under src/: tests/cpp unchanged and passing.
- The full tinytest suite on the shipped build, in one process, counted file by file: no failure, no file
  stopping, and at least the base build's count plus the new files' assertions less those removed with a
  law; the landing note gives the figures.
- On a reference build (`--preclean --configure-args=--enable-reference-build`): the four
  `test-reproducibility-*.R` files pass unchanged, and the three compares are bitwise, every scenario
  reporting identical draws, counted scenario by scenario with no `max |z|` line: 55 against the gaussian
  baseline in force, the BCF baseline in force, 11 against `multinomial-equivalence-80b1c8d4.rds`. In
  push 1 the BCF compare is first run against the old baseline, 14 identical and `glue_toggle` not, and
  then against the re-recorded one, 15. `glue_toggle` fails the statistical mode against the old
  baseline, as a changed posterior must (111 of 209 summaries beyond z = 3 on the stand-in); what shows
  its new values right is step 1.5's identity. In push 2 the compares are run after the first commit
  alone, 55, 15 and 11, and again at the push's end, 55, 15 and then 19 against the file with the four
  scenarios added.
- Every gate [exact-gates.yaml](../../.github/workflows/exact-gates.yaml) lists, in `quick`, each on its own
  exit status; in pushes 1 and 2 the four BCF scripts and forest-law-exact.R also in full.
- `inst/include/dbarts/dbarts.h` has no diff and `tools/check-api-hash.sh` passes.
- The pair script (below), old side on the base build, new side on the push's.
- `Rscript -e 'lintr::lint_package()'`, `air format --check .`, `Rscript tools/check-rc-codoc.R .`,
  `Rscript tools/check-win-drift.R .`, `Rscript tools/check-doc-freshness.R .`,
  `Rscript benchmarks/R/mutation-battery.R verify-anchors`, each on its own exit status;
  `pkgdown::check_pkgdown(".")`; `R CMD build` with every vignette rebuilt and `R CMD check --as-cran` on a
  tarball from a clean copy.
- The consumers, each suite whole against a private install of the push, none failing: bartCause on
  dbarts-1.0 (1412 expectations at its last run; with step 1.8's edit from push 1, and once without it on
  push 1, failing only where `amplitude.prior.variance` is passed), stan4bart on bartcore (582), treatSens
  on dbarts-1.0 (306), bairrtt on main (207), the last three unedited.
- A hot path is touched once, by push 2's first commit: `bench-sampler.R compare` on a quiet machine is
  the maintainer's to run, on that commit against its parent.

The pair script. Each row is one seeded model on both builds: sampler fits compare the train draws, sigma
and the coefficients, `bart` fits `yhat.train` and sigma. Rows 01 to 18 and 22 were run for the first
plan on the tip against the stand-in, which is the law at once: identical; rows 19 to 21 differed.

| | model | push 1 | push 2, first commit | push 2 |
|---|---|---|---|---|
| 01 to 03 | two forests, default; an sd on both; three with a factor third | identical | identical | identical |
| 04 to 09 | probit; logistic; weights and an offset; a term under `subset`; `bart` with two chains; a logical | identical | identical | identical |
| 10 to 14 | a swap to another factor; a write; a copy; a reload; `setResponse` | identical | identical | identical |
| 15 to 18 | no forest without a basis; forest weights and a mask; four forests; counts and a tree prior on a factor forest | identical | identical | identical |
| 19 | a numeric column, default and stated; a 0/1 number; a factor beside a number | identical | identical | differ |
| 20, 21 | a held factor; a held forest with no basis, default and stated | differ | as push 1 | as push 1 |
| 22 | a held forest with no basis, default, probit | identical | identical | identical |
| 23, 24 | against the base build's own spelling of the same prior: a held two-level factor, second forest, at s, the base build stating 0.674 s; a held forest with no basis at one unit, the base build stating any sd | equal to 1e-10 | | |
| 25 | a factor with `amplitude.prior.variance = 2` | refused | refused | refused |
| 26 | a held factor third; a held forest with no basis second | created (refused on the base build) | | |
| 27 | a held numeric column | refused | refused | created |
| 28 | dose in units that put its values near 1e-9, default | identical | NOT identical: the forest sees its rows | the fit of dose in any other units, to the bit |

Rows 23 and 24 are two changed sequences against a prior the base build can state: run for the first
plan on the tip and the stand-in, four such pairs under gaussian, probit and logistic were identical.
Where the base build has no such spelling (row 19) the oracle is the exact gate and the tests/cpp
literals.

Mutations, each expected to fail the named test and no gate before it. Apply, install with `--preclean`
where it is the engine's, run, record the failing count, revert, `touch` the file.

- push 1, the kind not carried from the data door (every basis numeric there): step 1.2's kind test, and
  step 1.3's held factor through `bases`;
- push 1, the kind lost on a copy or a reload: step 1.2's second half;
- push 1, a held value left by position: step 1.1's six orders, and step 1.3's third forest;
- push 1, a held forest with no basis second held at 0: step 1.3's plain forest written second;
- push 1, the 0.674 left on a held factor: step 1.1's leaf scale, step 1.3's `k.scale`, and
  forest-law-exact.R (0.10 and 0.21 at its two stated arms on the draft);
- push 1, a held forest with no basis left at the response's scale: the same three (0.033);
- push 1, its default taken as 1 under gaussian: step 1.3's default;
- push 1, a leaf's sd not divided by the root of the tree count; the default taken for K - 1 forests:
  the gate's arm of three trees (167 and 70 on the draft);
- push 1, `amplitude.prior.variance` taken and ignored: step 1.3's unused-argument pin;
- push 1, a write to a held forest with no basis sent to the half-Cauchy median: step 1.1's write;
- push 1, a held block installed from a state: step 1.1's state;
- push 1, one of forest-defaults-by-kind's interim refusals left standing: step 1.3's creations;
- push 2, the bound left absolute: step 2.1's twin, and step 2.8's units twin downward;
- push 2, the bound taken relative to the largest multiplier of every forest, or of the last sweep: step
  2.1's row at 2^-30 of its own forest's largest;
- push 2, a zero multiplier kept in the update: step 2.1's held block;
- push 2, the default per unit (the scale read as 1): step 2.5's `sqrt(2 / K) * U / sd(dose[kept])`, and
  the gate (0.06 to 0.85);
- push 2, a scale 3 percent off: the gate, every arm (4.6 to 10.6 on the draft);
- push 2, n for n - 1: step 2.2's 0.25 and 8, and step 2.5's 1e-12 (the gate sees it on the drawn arms
  only just);
- push 2, the scale over the nonzero rows; over the rows of positive weight; before `subset`; with
  weights: step 2.2's half-zero column, step 2.5's "which rows", and the gate's weights arms;
- push 2, one pass for the mean: step 2.2's 1e7 + i;
- push 2, "constant" decided by the computed sd: step 2.2's 0.1, 61.7 and 1e8 + 0.1, and step 2.6's
  cleared record; the near-constant bound dropped: 61.7 + 1e-9 sin(i);
- push 2, a numeric coefficient's variance left at 1/2; the row norm left in; the column centred: step
  2.3's literals and its closed form, and the gate (0.11; 0.51);
- push 2, a stated sd divided by the scale: step 2.3's literals, and the gate's stated arms (0.14 to 0.65);
- push 2, the scale taken again at a re-creation; at a swap: step 2.5's "once";
- push 2, the record kept where a declaration replaces the basis: step 2.5's "other data";
- push 2, the record written back on every creation: the caller's object in step 2.5;
- push 2, an unread record checked; a read one not: step 2.5's "carried and wrong";
- push 2, the record dropped by a write: step 2.5's "stated, then put back";
- push 2, a numeric coefficient started at 0 as the second forest: step 2.3's literals;
- push 2, a held numeric column at 0: step 2.6, and the gate's held arms (0.27);
- push 2, the reader's `sd` taken from the statement with the engine's unread: step 2.7's "through the
  engine";
- push 2, a numeric forest's entry left unnamed; `extract` returning a list of one: step 2.7;
- push 2, the constant-column check made before `subset`; the prior-predictive refusal dropped: step 2.6;
- push 3, a column's printed sd read from the record; the per-sd figure from the column in force where
  the default came from another; the range left out; a label by position; a mean not passed through
  `zapsmall`: step 3.2.

## NEWS

No new item: forests, `forest()`, its `sd` and `amplitude`, `$setForestBasis` and every reader here are new
in 1.0-0 and nothing released changes. The "Multi-forest models" item is reread and changed only where it
names `amplitude.prior.variance` or the 0.674.

## What this leaves for later, and why each is an addition

Each requirement with what holds it; "test" is the step above that fails if it is broken.

| later | what it needs of 1.0 | held by |
|---|---|---|
| `$setForestBasis(updateBasisScale = TRUE)` (dec-B265, dec-B276) | the argument does not exist, so no call states it | R's unused-argument error, pinned in step 2.6 |
| | with no argument a swap keeps the scale: `FALSE` is today's meaning | step 2.5's "once" |
| | a record a call can rewrite, and the arithmetic to restate a forest from a scale | the slot of step 2.5; `basisColumnScale` and `applyForestLaw` of 2.2 and 2.3 |
| | no state, no shipped signature: the update is a virtual of the engine's own facade and one more bridge argument | Constraints |
| `leaf.prior = normal(sd = )` on `forest()`, and `normal()` for the default back (dec-B266, dec-B276) | `leaf.prior` on `forest()` is R's unused-argument error | written-surface's pin |
| | the reader already says "stated" or "the default" | forest-sd-unit's `sd.stated`; no test here pins the class of `leaf.prior` |
| | a stated forest keeps its record, so a default put back has the creation's scale | step 2.5's "stated, then put back" |
| a heavy-tailed coefficient, `amplitude = student(3)` (dec-B262, dec-B271) | `amplitude` takes `fixed()` and nothing else | written-surface's pins |
| | one row of `forestLaw`: a drawn numeric or factor block under the forest with no basis's law | the law is one function |
| | the state already carries one variance a forest and installs it where a size is half-Cauchy | Context; no edit |
| | the reader has an `amplitude` entry to say it in | step 1.3 |
| a forest for each level of a factor of three or more (dec-B281, TODO `factor-basis-per-level`) | the call is refused today | the kind-by-class slice's pin |
| | each such forest is a two-level one, which every reader and the block already describe by label | nothing more of this slice |
| a basis of several numeric columns, either meaning (dec-B282, TODO `several-column-basis-meaning`), and a size for each column (dec-B272) | the call is refused today, and so is an `sd` of length two | the kind-by-class slice's pins; written-surface's |
| | the reader's `sd` and `basis.scale`, and `extract`'s entries, are vectors by column, named by the column | step 2.7 pins the name on the one column |
| | `data@basis.scale` holds a vector for each forest | step 2.5 |
| | the engine is told a forest's kind and never reads it off a width | step 1.2; "The engine and 1.0-0" |
| | no stored state and no compiled interface holds a forest's sd or scale | Constraints |
| | not built now, each behind the bridge when it comes: a prior scale for each coordinate of a block, a default for each column, a law and a held value for a numeric block of several columns | "The engine and 1.0-0" |

## What waits on what

- On push 3 of written-surface: a basis's column name, which names a numeric forest's `sd`; the labels
  that name the reader's and `extract`'s lists and the texts; every call here is in its spelling.
- On forest-defaults-by-kind: a forest with no basis may stand anywhere, so step 1.3 holds one written
  second; its interim refusals of held shapes, two of which push 1 lifts; selection by label in the
  reader; its line of tree counts, which the block follows.
- On forest-sd-unit: the unit, the statement the engine is handed and the one function that holds the
  law, the four-number record, the reader's `leaf.prior`, `sd` and `sd.stated`, `extract`'s list by label
  and its entry for one forest, a write that states. This slice changes rows of that function, adds the
  kind to the statement and names one entry.
- On forest-kind-by-class: the record of a factor's levels and the one function that gives a kind; its
  second push, after which every forest is one numeric column or two levels and none changes width, so
  that what moves in push 2 here is a number that was meant as one; bartCause's factor line.
- To recheck once each has landed: the stand-in is of no further use for counts, most of what it failed
  being repaired or removed by then; run the suite on the push's first build and count. The pair
  script's rows. That no baseline scenario has gained a number or a hold. That `samplePriorPredictive`
  still cannot predict a sampler of several forests (if it can, it returns draws with the coefficients
  at their starting values, and step 2.6's refusal moves to push 1). That the suite still holds only the
  two held shapes forest-defaults-by-kind kept. The names of the functions cited here.
- With [leaf-conversions.md](leaf-conversions.md): it edits the state's reader and writer, `Chain`'s
  install and `setData` paths and the leaf models; this slice edits the K-forest constructor, the forest
  setters, the amplitude combiner's constructor and one of its methods, and the bridge's creation path,
  and no line of the state. Serial for the shared files, in either order, each with `--preclean`.
  cross-family-state-install has landed, and the removal of its family record (dec-B283, TODO
  `state-family-record-removal`) edits the bridge's state code and nothing of this slice's.
- With the state-frame arc's data record: `basis.scale` is a fourth record of its kind, data derived once
  and then held, and follows its rule (an object that carries one: used). It stays on the data object
  (dec-B265, confirmed by the maintainer). This slice adds the slot and its writer itself and waits for
  nothing; whichever lands second folds the two read-backs into one. Two things that arc must not do to
  it: write it after a state install, which never moves it, or take it off the object a copy hands to
  creation, which reads it there; step 2.5's copy after a swap is the guard, and step 2.9's scenario
  with a swap, a write and a copy the recorded one.
- With the control-migration arc: its first slice moves the forests' record to the model, whole, as
  forest-sd-unit leaves it. Serial; this slice first keeps the posterior-changing gates away from slices
  gated bit for bit.

## Out of scope, and where it goes

- `updateBasisScale`, the long form on `forest()`, the heavy-tailed law: after the merge (TODO
  `forest-after-merge`, `coefficient-law-opt-in`). A forest for each level, and several numeric columns
  with a size for each: after 1.0 (dec-B281, dec-B282, dec-B272).
- `setData` on a sampler of several forests: refused, as today. The rule for the day it opens goes in the
  design note: a new object's record is used, what it lacks is taken.
- Holding or releasing a coefficient on a live sampler: no setter exists; `setModel`'s refusal of a model
  that differs in it waits for the forests' record to be on the model.
- In TODO already, added with this amendment: a state from a sampler that draws a two-level factor's
  coefficients installs into one that holds them, in silence.
- To TODO as new entries: a drawn coefficient starts by its forest's position (1, and 0 for the second
  forest's first), a start by kind being a change of draws for a model of three forests;
  `samplePriorPredictive` for several forests, which needs the coefficients drawn from their prior; a
  consumer that creates a sampler of several forests through the shipped header gets a scale taken at
  each creation, there being no R object to record it on (none does).

## Ruled after this plan was amended

To be worked into the steps when the plan is rechecked before push 1 is built; where a step above says
otherwise, this section stands.

- dec-B290. Under a gaussian response a held forest with no basis that states no sd takes the default
  `bart()` gives the one forest of an ordinary fit on the same response, a quarter of the response's
  range over the rows kept at k = 2, about 1.4 sd(y) for a normal sample of 300, and not the 2 sd(y) of
  the steps above. The drawn forest with no basis keeps 2 sd(y) as the median of its scale. The engine
  holds the rule, as it holds every default, and takes it from the code that gives `bart()` its own, so
  the two cannot drift; the help states it as "the size bart() gives its forest". The default under
  probit and logistic is not ruled and stays as planned. bartCause's `bcf(update.a = FALSE)` states
  `sd = 2 sd(y)` on its first forest the day push 1 lands, bcf's own convention (`sd_control`), so its
  fits keep that prior: one line in R/bcf.R and its help.

      forest(x1 + x2, amplitude = fixed())            # own sd range(y) / 4, as bart()'s forest
      forest(x1 + x2, amplitude = fixed(), sd = 2 * sd(y))   # what bartCause's bcf passes

- dec-B291. The default sd of a forest with a basis does not depend on the number of forests: it is 1 on
  the response's scale (sd(y) under gaussian, the link's unit otherwise), per standard deviation of the
  column for a number, where the steps above have sqrt(2 / K). At two forests the number is unchanged.
  Every step, test, gate arm, help sentence and refusal text that carries sqrt(2 / K) is reworked; the
  same holds for [forest-sd-unit.md](forest-sd-unit.md) and
  [forest-defaults-by-kind.md](forest-defaults-by-kind.md) wherever they state the default.

## Calls made in planning

The coordinator's calls on the critique, each to be recorded in the ledger at landing:

- The row bound becomes relative to the forest's own largest absolute multiplier in the sweep (the
  critique's finding 1). It is push 2's first commit, gated bit for bit alone, with the bench and the
  three baseline compares; the help promises a fit that does not depend on a column's units only with
  it in.
- Every held shape that is wrong on the tip is refused in the interim, by width and position, in
  forest-defaults-by-kind's first push (finding 2); push 1 here lifts what remains.
- The gate decides by z, with a floor of 0.0005 and 32 seeds of 25000 (finding 3, with the previous
  planner's correction): the first plan's fixed 0.01 passed a scale 3 to 10 percent off. One arm has
  three trees a forest and K = 3. The three gaussian BCF gates leave the list of oracles; they cannot be
  made to see the held law. The held arms and the numeric arms live in one script. One prior-size test a
  family runs always.
- Finding 4, two refusals that told the user to state one sd for columns in different units: both
  concerned several columns and are gone with them. What is kept of it: the constant column's text says
  what to state, and the help says what a stated sd is per unit of.
- Four scenarios enter the BCF baseline when push 2 records (finding 5), the critique's two-column one
  dropped.
- "Constant" is decided by comparing values, in R and in the engine, and a column whose standard
  deviation is under about 1e-8 of its largest absolute value is refused with the constant's text
  (finding 6); the bound is written 2^-26, the engine's own constant for "the same number".
- The first push as first planned is folded into forest-sd-unit, and the kind's route to the engine
  moves to the held push, the first that reads it (finding 7).
- Of the critique's questions: (a) a held forest with no basis keeps the planned default and is listed
  above as a question still to be put; (b) which shapes may be held is settled by dec-B281 and dec-B282;
  (c) `extract` for one forest returns the entry itself, as
  [`selectForests`](../../R/generics.R) does for every other quantity, and a list for several; (d) a
  factor's entry is one unnamed number, the size of a difference between levels, which belongs to no
  level; (e) the reader gives `forest()` for a default with the number beside it, and tests pin `sd`,
  `sd.stated` and the round trip, not the entry's class; (f) a constant beside a column is refused by
  dec-B282; (g) is not a question: the scale stays on the data object (dec-B265, confirmed by the
  maintainer; the state-frame design classes it as a data record); (h) the printed block stays before
  the merge, and prints the column's range beside its mean and sd.
- Cut: the refusals and tests for objects only a development build wrote (the bridge's check of the
  record's length stands), and a fit's stored mean and sd for print. The scale's computation stays in
  the engine, which holds every default and is never R-specific.

The planner's own, on the amendment:

- Nothing is built for several columns. dec-B272 has the engine hold one scale a column from this slice;
  under dec-B282 no 1.0 model has a second column, and what that sentence protects, no stored state and
  no compiled interface changing later, holds by the state and the header holding no sd ("The engine and
  1.0-0"). Kept instead, at no cost: the shapes a 1.0 user sees (a vector by column in the reader, in
  `extract` and on the data object) and the kind stated to the engine.
- The kind is stated to the engine though width alone would tell the three 1.0 shapes apart. A law read
  off a block's width would have to be unpicked the day a numeric block has two columns, and would make
  the engine's reading of a block depend on what R happens to refuse.
- The engine is not narrowed to what R accepts: a drawn block of more than two levels keeps the law it
  has. Removing it would be work to undo nothing R can reach, and tests/cpp would lose fixtures that are
  not this slice's.
- Three pushes, cut between the held coefficients, the numeric law and the printed block. The first cut
  lets the held law land on the gates that exist and gives the numeric law a hold that already goes by
  kind; its cost is that a held numeric column stays refused for one push more.
- The bound's change is one commit inside push 2 and not a push: it is six lines, it is what the push's
  twin tests stand on, and a separate battery for it would repeat the three compares the push runs
  anyway. It is gated alone all the same, because it is the one change here that could move a draw of a
  model with no number in it.
- A lone z between 4 and 5 reruns its arm once on other seeds. At about a hundred functionals a run a
  bound of 4 trips by chance in a few runs of a hundred, and a gate that fails then teaches people to
  rerun it without reading it.
- The gate reads nothing off the sampler but the unit and the shift, and takes the scale from R's `sd()`.
  The alternative, reading the scale back, would pass any wrong scale.
- A scale is taken only for a forest that states no sd at creation, and a later write leaves the record.
  Kept, a default put back later has the creation's scale. No width changes, so no entry goes stale.
- A record is dropped where a declaration replaces a forest's basis on a data object. A hand-edited object
  keeps its record, as the state-frame design rules for its three: to take a new one, clear the entry.
- The mean is taken first and the deviations squared, as R's `sd()` does: on a column far from zero one
  pass differs from R, and the printed line would name a scale the user cannot reproduce. The first
  plan's test of it, 1e8 + 1e-8 k, is a column the near-constant rule now refuses; 1e7 + i replaces it.
- The scale is unweighted and counts rows of weight zero, as the response's scale does: one estimator for
  both, and the one R's `sd()` gives for the column the caller sees.
- A numeric forest's `sd` is a named vector of length one, and a factor's an unnamed number. dec-B275
  rules a vector by column, and dec-B282 keeps that shape "with one number in each".
- The reader says "factor" where the kind-by-class slice's function says "levels": one is the word the
  user wrote, the other names the columns.
- Gone from the reader: the four entries of the tip's channels; the divisor is in `prior.sd.of`'s "(prior
  median)".
- A state holds nothing of the law, so none is refused and none is re-read. A state stored on a build that
  held a forest's coefficients at other values installs its trees under the recipient's, by the standing
  rule; a check there would refuse what
  ["fixed amplitudes"](../../inst/tinytest/test-state-not-model.R) pins as taken. It is in TODO.
- Drawn starts stay where they are for a forest with no basis or a factor. A start by kind would be the
  cleaner rule and would move the draws of every model of three forests.
- A held numeric coefficient with no sd takes the drawn default, as the design has it: right at ordinary
  slopes, short under probit when the slope is two or three times larger, which the help says.
- The print's "per sd" figure stays though it is one number for every default line: it is what shows a
  stated sd left unchanged when its column was rescaled (dec-B266). The range is what shows one wild
  value.
- The verbose summary gains the block from R; the engine's own summary, which prints the first forest's
  model, is left to the TODO entry forest-sd-unit adds.
- Measured for this plan on stand-ins, not on the slice: the law's rows, the 20 pairs, the suite's counts
  and the BCF compare come from the design's prototype ported to the tip, which holds the kind as an
  attribute, the sd in units of L and six numbers a forest; the gate's figures from a draft run against
  the prototype's own build; the relative bound from the critic's six-line copy of the stand-in. Read
  off push 3's build for this amendment: the held values by position, the counts of held forests and of
  constant columns, the reader's and `extract`'s shapes. Read and not run: the three baseline compares
  and the bench on the relative bound, the full gaussian compare, the snapshot files, sanitizers, the
  consumers' suites. Nothing of the three pushes is built.
