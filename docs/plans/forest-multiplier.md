# forest-multiplier: what a forest's coefficient is, held or drawn, and what a numeric multiplier's sd is the size of

Status: PUSH 1 LANDED 2026-10-17 (f9bc3636..5a62933e, dec-A216; [Push 1 as built](#push-1-as-built), [Landing note](#landing-note)); PUSH 2 LANDED 2026-10-17 (a51ac7fe..5a1e9dca, dec-A219; [Push 2 as built](#push-2-as-built), [Landing note, push 2](#landing-note-push-2)); push 3 PLANNED,
rechecked 2026-10-11 against 37b59ea0, the tip of the sd unit slice
([Recheck against bartcore](#recheck-against-bartcore)); its two open calls settled and its cites
brought to forest-kind-by-class as built on 2026-10-17, with push 1 gaining the unit on a constant
response. Rulings: dec-B253 and dec-B246 as revised;
dec-B260 to dec-B265; dec-B269 as dec-B280 and dec-B282 revise it; dec-B271, dec-B272, dec-B275,
dec-B276, dec-B281, dec-B282; dec-B290, dec-B291 and dec-B293 for the defaults; dec-B303 with dec-B302
for the rows a column's scale reads; dec-B307 for a constant column; dec-B346 and dec-B452 for a swap;
dec-B353 and dec-B354 for a held coefficient; dec-B362 and dec-A207 for a re-derived response scale;
dec-B376 for what a fit records of a hold; dec-A171. Under dec-B354 a held coefficient is 1 and the
basis carries any coding: `fixed()` on a forest refuses a value other than 1, and the engine no longer
holds by forest position. The two calls the recheck left open are settled, each as the option the steps
were written to ([Settled calls](#settled-calls)). Follows [forest-defaults-by-kind.md](forest-defaults-by-kind.md),
[written-surface.md](written-surface.md) and [forest-sd-unit.md](forest-sd-unit.md), landed; and
[forest-kind-by-class.md](forest-kind-by-class.md), whose push 1 landed and whose push 2 is built
(32418b6e) and in its landing gates. Amended 2026-10-07 after its blind critique (verdict: build after corrections) and after dec-B281
and dec-B282, under which a basis is one numeric column or a factor of two levels: three pushes where
there were four, and nothing for several columns. The counts of tests, lines and blocks below are those
of 2026-10-07 and are retraced on the first build after forest-kind-by-class's second push, which
rewrites the tests they count.

One plan, three pushes. The engine's change of rule and the R surface that states it are not planned
apart: at every tip the help must say what the engine does and the reader must report the number the
engine uses, and a plan for the engine alone would leave a tip that does neither.

agent: pushes 1 and 2 (the rule): opus implementer for the engine, the bridge and the R code, sonnet for
the respelled tests, the benchmark scripts and the help once the code is fixed, opus reviewer told to
refute. Push 3 (the printed block): sonnet implementer, opus reviewer. The reason for opus on 1 and 2:
every slip is a prior off by a factor and no message. A column's scale is taken at two places of the
engine today, the K-forest constructor of [`Chain`](../../src/bartcore/chain.hpp) and
[`Chain::setForestBasis`](../../src/bartcore/chain.hpp), and every creation from R ends at the first,
by the nine routes this plan first counted as paths: six creation doors (seven as
forest-kind-by-class counts them, a `forests` list on a formula once for a basis written as code and
once for a value), a copy, a reload and a sampler made again from its own control, model and data. On
any of them the scale can be taken from the wrong rows, or taken again where the record should be handed
over, and a held value can go by position on one route and by kind on another; a stand-in for the
slice, built with care, had three such faults (Context).
rng: stated per push and per call sequence, as [RNG classes and their gates](README.md#rng-classes-and-their-gates)
defines the classes.
- Push 1, its first commit, the default of a forest with a basis (dec-B291): POSTERIOR-CHANGING for a
  model of three or more forests in which a forest with a basis states no sd, a factor or a number,
  drawn or held: its default is one unit where it was sqrt(2 / K) of one. Bit for bit at two forests,
  and wherever every forest with a basis states its sd. The oracle is an identity: the base build's
  fit of the same model with each such forest stating the unit (1 under probit; under gaussian and
  logistic the `sd` the reader gives a default factor forest of a two-forest sampler on the same
  response), identical draws.
- Push 1, the holds: POSTERIOR-CHANGING for a model that holds a coefficient, in two sequences: a held
  two-level factor (the forest's own sd is s where it was s / 0.674); a held forest with no basis, under
  every family and whether or not it states an sd (its own sd is s, and under gaussian a quarter of the
  response's range where none is stated, dec-B290, where it was the response's scale whatever s; under
  probit and logistic at its default of one unit the fit keeps its bits, dec-B293).
  NEUTRAL for every model whose coefficients are all drawn, a numeric multiplier included, the first
  commit's sequence apart. Accepted
  where refused, forest-defaults-by-kind's two interim refusals being lifted: a held forest with no basis
  that is the second forest, at 1; a held two-level factor that is not the second, at 0 and 1. Refused
  where accepted: `amplitude.prior.variance`.
- Push 1, the unit on a constant response (step 1.9): POSTERIOR-CHANGING for a gaussian model of
  several forests whose response has one value over the rows its scale is read from: every forest's
  default sd was 0 and is now its default times a quarter of the response model's window of width
  one. Accepted where refused: an sd stated on such a model, and its `copy()`, reload and
  `new("dbartsSampler", control, model, data)`, which the record's anchor of 0 refused.
- Push 2, its first commit, the bound under which a row leaves a forest's update: bit for bit for every
  recorded scenario, every snapshot and every seeded pin of the suite, and gated so, alone. In law it
  changes a sweep in which some row's multiplier is under 2^-26 in absolute value and not under 2^-26
  of the forest's largest, or the reverse; that is every sweep of a numeric column whose values are
  below about 1e-8, to which the forest was blind, and otherwise a sweep where a coefficient has come
  within nine orders of magnitude of zero or of another level's.
- Push 2, the rule: POSTERIOR-CHANGING for every model with a numeric multiplier, drawn, at the default
  sd or a stated one. NEUTRAL for every model whose forests have no basis or a factor, held or drawn,
  stated or not. Accepted where refused: a held coefficient on a numeric column (dec-B353's refusal is
  lifted). Refused where accepted: a numeric basis that is constant, or constant to rounding, over the
  rows its scale is read from, on a forest that states no sd (dec-B303, dec-B307); with an sd stated it
  is created, as today.
- Push 3: no draw moves.
- Recorded scenarios. Push 1 moves one, `glue_toggle` of the BCF equivalence baseline, which holds the
  treatment forest's coefficients, and re-records it. Push 2 moves none and adds four to that baseline,
  so that the rule it builds occurs in a recorded draw. No other scenario of the three baselines holds a
  coefficient or multiplies a number, every recorded model of several forests has two, so the first
  commit moves none, and the four snapshot files fit no model of several forests
  (searched 2026-10-07, and for the count of forests by the sd unit's recheck at d5b637ee; searched
  again before the build, forest-kind-by-class rewriting one gaussian scenario in between; the BCF
  compare was run on a stand-in: Context).
window: pre-release, after forest-kind-by-class's second push and before the control-migration arc moves
the forests' record to the model. Serial with any other work in
[`ForestSpec`](../../src/bartcore/combiner.hpp), [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp),
the K-forest constructor of [`Chain`](../../src/bartcore/chain.hpp),
[`Chain::setForestBasis`](../../src/bartcore/chain.hpp), [`Chain::setForestMapSd`](../../src/bartcore/chain.hpp),
[`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp), [`parseData`](../../src/R_interface_bartcore.cpp),
[`bartcore_getLeafPrior`](../../src/R_interface_bartcore.cpp), [`forestParams`](../../R/model.R), the
multi-forest block of [`resolveSamplerSpec`](../../R/spec.R), [`reportLeafPrior`](../../R/dbarts.R) or the
sampler's `initialize`. bartCause's edit lands the day push 1 does.
budget: as rechecked 2026-10-11, ~5800 lines planned and a forecast of 9000 to 12500: push 1 ~2300
(forecast 3600 to 5000), push 2 ~2950 (4600 to 6400), push 3 ~520 (800 to 1100). The forecast is the
plan times 1.55 to 2.15: the sd unit slice was built at 1.56 times its plan (3646 lines against 2340),
and the four slices before it at 1.8 to 2.3. Against the figures below, push 1 gains the default of
dec-B291 (~180: one line of the engine, the pins of sqrt(2 / K) at three and four forests, the help),
the held default of dec-B290 (~130) and nine more test lines that state `amplitude.prior.variance`,
and loses about 50 where forest-kind-by-class's function and loop are reused; push 2 gains the rows of
dec-B303 with their fallback (~150) and loses a little where the tests' creations on a constant column
have gone from 55 to 31 (forest-kind-by-class's trace). Stops: push 1 at ~5000 or at ~900 of package
code, push 2 at ~6400 or at ~1100 of package code; past either, the push is cut at its last gated commit and what remains is
planned again. As amended 2026-10-07, and not recounted since:
~5400 lines changed, upper figure 8900. By push: 1 ~2050 (upper 3380), 2 ~2830 (upper 4670),
3 ~520 (upper 850). By layer, with the upper figure: engine 370 (610); bridge 230 (380); tests/cpp 600
(990); R 740 (1220); tinytest 1720 (2840); benchmarks 800 (1320); help 430 (710); design note,
architecture, TODO and the two indexes 510 (830). As first planned the slice was 7800 lines in four
pushes. Gone from it: the neutral first push, folded into forest-sd-unit (about 950, the kind's route
staying here); everything for several columns and for a width that changes (about 1300); the refusals
of development-build objects and a second gate script (about 300). Added: the relative bound with its
twins (about 150), four scenarios (80), an arm of three trees (60), the tests that run always (40).
Past estimates in this arc ran 1.5 to 2 times low; the upper figure is 1.65 times the estimate.

## Goal

Each kind of forest has one rule, and it is the rule the help states. A forest's `sd` is the size of what it
contributes for one step of its multiplier: exactly where its coefficient is held, as the prior median where
it is drawn. A numeric column is used as it was given, never centered and never divided; its coefficient is
normal; a stated sd is per unit of the column; and with no sd stated the default is one unit of the
response's scale for one standard deviation of the column, whatever the number of forests (dec-B291),
that standard deviation taken once, over the rows in the likelihood when the sampler is created
(dec-B303), and kept with the data. A held coefficient is 1 and the multiplier is what the basis says,
wherever the forest stands: 1 with no basis, 0 for the first level of a two-level factor and 1 for the
second, the column as given for a number (dec-B354). A held forest with no basis that states no sd has
the size `bart()` gives its one forest (dec-B290). `amplitude.prior.variance` is not an argument. A
default fit is the same fit whatever units a
column is in, down to units in which its values are of order 1e-12. The reader, `extract` and `print`
give a numeric forest's number under its column's name, and it is the number the engine draws under.

## Context

Measured at d02fe72c (written-surface pushes 1 and 2 landed; push 3, forest-defaults-by-kind,
forest-sd-unit and forest-kind-by-class not) on the shipped build, R 4.6.1, unless a line says otherwise.
Read again on 2026-10-11 at 37b59ea0, the sd unit's tip as landed (read at its re-gate tip, the same code), in the code and by four scripts run on its shipped
build (scratch/fmr/ in the recheck's worktree): a bullet that begins "At 37b59ea0" is that reading, and where it stands beside an older figure
the older one is history. Every count of tests, files and lines is of the older tip and is retraced
before the build. The spellings are the tip's but for a tilde on a basis
written as code, which the landed grammar no longer writes. L is the response's scale as the engine holds it: sd(y - offset) under gaussian, 1 under
probit, pi / sqrt(3) under logistic, the first read over the rows in the likelihood at creation
(response-scale-rows, landed); the sd unit slice makes it the unit a stated sd is divided by, so that s
below is a number of response units over L, and changes nothing else below. K is the number of forests. "Size" is the prior median of the absolute multiplier
times the prior standard deviation of the forest's own fit, from the engine's prior draws under a flat
likelihood (4000 sweeps a row, Monte Carlo error about 0.02); a held multiplier is exact. Fixture: 300
rows; x1, x2; a 0/1 z; w with standard deviation 3.15; an age with mean 50. Rows for a factor of three
or eight levels and for two numeric columns were measured too; both shapes are refused before this slice
(dec-B281, dec-B282) and are left out.

- The suite. At the tip 232 files; on push 3's build (4eeaf03d, traced 2026-10-07) 233, of which 230
  without the three that ask for more than two threads give 17121 results, none failing. At 37b59ea0
  there are 257 files; not run for this recheck.
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
- At 37b59ea0 the table stands, row for row, with s in response units over L. The engine resolves it in
  one function, [`forestMapEntry`](../../src/bartcore/combiner.hpp), from what each forest states
  ([`ForestStatement`](../../src/bartcore/combiner.hpp)), and
  [`defaultForestSd`](../../src/bartcore/combiner.hpp) still gives sqrt(2 / K): the sd unit slice left
  dec-B291 to this one. Run: three forests on a response of sd 1.124, a forest with no basis, dose and a
  two-level factor, give `sd` 2.247, 0.917 and 0.917 (2 and sqrt(2 / 3) of the unit); the numeric
  forest's `k.scale` is 68.17, its sd over 0.674 and over the median absolute dose, 0.01997; a held
  forest with no basis reports `sd` 2.247 and `k.scale` 1.124, the response's scale whatever s; a
  two-level factor held with `fixed(1)` reports `k.scale` 1.667, its sd over 0.674.
- `amplitude.prior.variance = v`, on a forest with a basis: the coefficients are N(0, v) and the forest's own
  sd does not move, so the size is sqrt(2 v) times the table's: 1.95 at v = 2 and 0.52 at v = 1/8 on a factor
  (2.04 with `sd = 0.7` beside it), 1.03 / c at v = 1 on a number. Nothing on a held forest. Refused on a
  forest with no basis. Stated on 31 lines of 6 test files, in 5 benchmark scripts and by bartCause.
  At 37b59ea0: an argument of `forest()` still, handed to the engine as written
  ([`ForestStatement`](../../src/bartcore/combiner.hpp)); on 40 lines of 7 test files, in the same 5
  scripts, in R/model.R and R/dbarts.R, and on two help pages.
- A held value goes by position. [`AmplitudeState`](../../src/bartcore/combiner.hpp) starts at 1 for the
  first forest and (0, 1) for the second, and [`rebuildAmplitudeLayout`](../../src/bartcore/combiner.hpp)
  carries those by position and fills every other coefficient with 1. Read off the engine, and again off
  push 3's build: a two-level factor held as the third forest, or as the first where no forest is
  without a basis, is (1, 1): no contrast, a function added on every row (size of the difference 0.000);
  a forest with no basis held as the SECOND forest is held at 0. A drawn coefficient starts at the same
  values. forest-defaults-by-kind refuses those two held shapes, by width and position, until this
  slice's first push; the suite holds 49 forests on push 3's build, every one a forest with no basis
  first or two columns second.
- At 37b59ea0 a held value still goes by position, in three places of the engine: the starting vector
  of [`AmplitudeState`](../../src/bartcore/combiner.hpp) (1; 0 and 1),
  [`rebuildAmplitudeLayout`](../../src/bartcore/combiner.hpp), which carries by position and fills
  with 1, and the two-forest reading of [`restoreGlue`](../../src/bartcore/combiner.hpp), which passes
  by a held block. R keeps the wrong shapes out by width and position:
  retired: [`refuseHeldShape`](../../R/model.R) at creation (one numeric column, dec-B353; a forest with no
  basis written second in a `forests` list; two columns that are not second; three or more) and, at
  a swap, retired: [`refuseHeldWidthChange`](../../R/model.R), which forest-kind-by-class's second
  push removed with the last swap that could change a width. Run: of three forests, a forest with no basis,
  dose and a factor, all drawn, the coefficients start at 1, 0 and (1, 1), so a drawn numeric
  coefficient that stands second starts at 0; a formula puts the forest with no basis first wherever
  its term is written, so the refusal for the second place is met through a `forests` list only.
  `fixed()` with a value other than 1 is refused already
  ([`validateForestAmplitude`](../../R/model.R)), in words that point to `sd` and not to the basis; a
  fit records a held coefficient as 1 in `fit$fixed$amplitude`, named `forest<i>` (dec-B376). The
  suite's held forests by forest-kind-by-class's trace: 70, of which 38 two levels second, 27 with no
  basis and 5 two numeric columns second, the last refused by that slice.
- The row norm. c is taken over the rows the sampler holds, after `subset` (1.976 for rows 1 to 200, R's
  median of those rows), over the rows where the basis is not zero, unweighted; rows of weight zero and rows
  masked later count. It is taken again at every swap: `$setForestBasis(2, 3 * w)` divides the forest's own sd
  by 3 (2.493 to 0.831) and leaves the coefficient. `copy()`, a reload and
  `new("dbartsSampler", control, model, data)` take it from the column then in force, so after a swap to
  w / 10 all three have ten times the creation's sd. Nothing of it is recorded.
- At 37b59ea0 the row norm is as above and is taken at the two places the header names
  (retired: [`basisRowNorm`](../../src/bartcore/chain.hpp)), over every row, while the unit of an sd is read
  over the rows in the likelihood ([`GaussianResponse::scaleRows`](../../src/bartcore/model.hpp)). Run:
  with ten rows of weight zero whose dose is 5, the row norm is the median over all 150 rows (0.02092;
  over the rows in it would be 0.01957) and the unit is the sd of the 140 (1.2496; 1.2334 over all);
  a swap to 3 dose divides `k.scale` by 3 and a copy keeps the swapped figure; a column constant at
  61.7 with no sd stated is created, `k.scale` 0.0270; a sampler created with every weight zero warns
  that the response's scale was taken from all rows. A probit or logistic response keeps no record of
  which rows are in, there being nothing it reads from them.
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
  At 37b59ea0: `updateScale = TRUE` is refused on a sampler of several forests by `setResponse`,
  `setOffset`, `setWeights` and `setActiveRows` alike (dec-B362, dec-A207; the first and third run),
  `test =` beside a `forest()` term is refused (run), and `$setForestBasis` has no `updateBasisScale`
  (run: R's unused-argument error).
- The prior-only path. `sampleTreesFromPrior` with `sampleLeafParametersFromPrior` draws each forest's own
  fit at the sd the reader reports (three forests over 3000 draws: 0.997, 1.027, 0.987 of `k.scale`) and
  leaves every coefficient where it was. [`samplePriorPredictive`](../../R/dbarts.R) on a sampler of several
  forests stops inside `predict`, on push 3's build too, in a text about test columns. At
  37b59ea0 it still does (run: "this sampler's forest amplitudes have no off-sample basis"); dec-B361
  has it draw under the prior the sampler runs under where it draws at all.
- A state. Its glue block is K, the widths, every coefficient, held ones included, and one variance per
  forest ([`readAmplitudeGlue`](../../src/R_interface_bartcore.cpp)). An install takes a forest's
  coefficients where the recipient draws them and its variance where the recipient's size is half-Cauchy
  ([`restoreGlue`](../../src/bartcore/combiner.hpp)); a held block in the state is passed by, the rule
  ["fixed amplitudes"](../../inst/tinytest/test-state-not-model.R) pins. So a state from a drawn sampler goes
  into one that holds (0, 1), the trees arriving and the coefficients staying (0, 1), "exact" `TRUE`. A state
  holds no sd, no scale, no row norm and no kind. (The family record it had at d02fe72c is gone:
  dec-B283, landed.) Every install puts the chain's parameters in as stored, never converted
  (dec-B418), so a numeric coefficient from a sampler created on another column arrives as the number
  it was.
- The engine, as the sd unit slice built it (37b59ea0). [`ForestSpec`](../../src/bartcore/combiner.hpp)
  carries a [`ForestStatement`](../../src/bartcore/combiner.hpp): an sd in response units or none, a
  multiple of the anchor for fixtures, a hold and, where `amplitude.prior.variance` is stated, a
  variance. [`forestMapEntry`](../../src/bartcore/combiner.hpp) resolves the tip's rule from that, from
  whether the forest has a basis, from K, the family, the unit and the row norm, once for each
  construction and again at a write ([`Chain::setForestMapSd`](../../src/bartcore/chain.hpp)); the
  forests' record holds four numbers a forest. All the engine knows of a kind is that one bit, and
  retired: [`ForestSpec::createdWithoutBasis`](../../src/bartcore/combiner.hpp), gone with push 1, with `hasBasis` in the
  forests' record, keeps it for a forest created with no basis and given one since, a case dec-B452
  closes. The K-forest constructor of [`Chain`](../../src/bartcore/chain.hpp) still
  divides by retired: [`basisRowNorm`](../../src/bartcore/chain.hpp)
  ([`mapLeafScale`](../../src/bartcore/chain.hpp)). [`ForestAmplitudePrior`](../../src/bartcore/combiner.hpp)
  holds one variance for a forest's whole block, which is all one numeric column or two levels need.
  The bridge's swap still drops the engine's answer
  ([`bartcore_setForestBasis`](../../src/R_interface_bartcore.cpp)).
- The shipped header. [`DBARTS_C_API_LIST`](../../inst/include/dbarts/dbarts.h) creates no sampler from a
  forest's statement and reads no forest's spread, scale or hold; a per-draw callback sees the coefficients
  as `glue`. Nothing in this slice changes a signature there.
- The reader, as the sd unit slice built it (37b59ea0, run). `$getLeafPrior(f)` on a forest of several gives `leaf.prior`
  (the statement), `sd`, `sd.stated`, `leaf.model`, `prior.sd.of` ("amplitude scale" or "forest total"),
  `prior.mean`, `k.scale`, `response.scale`, `response.shift`,
  and `amplitude.prior.scale` or `amplitude.prior.variance`, `leaf.scale.factor`, `leaf.scale.divisor`,
  `basis.row.norm`; with no forest, a list named by label. The statement does not say that a
  coefficient is held. `extract(type = "leaf.prior.sd")` is the list of
  the `sd` entries by label, and the entry itself for one forest; `extract(type = "k")` is 1 for every
  forest, named `forest<i>`. `print` of a fit gives the tree counts by forest and nothing of any
  forest's sd, as `show` of a sampler and the verbose summary do; `summary()` of a fit has one line
  under its table, `(Fixed, not sampled: leaf.prior.sd[forest1] = 2.2471, leaf.prior.sd[dose] = 0.9174,
  leaf.prior.sd[zf] = 0.6000)`, the reader's `sd` entries by label.
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
  seeds and at -1.49 at 64 seeds of 50000: noise. The draft's oracle is the rule of 2026-10-07: it has
  sqrt(2 / K) in a default, 2 sd(y) for a held forest with no basis and a scale over every row, so
  that "the scale taken over the rows of positive weight" above is now the rule and "over every row"
  the fault (dec-B303). The oracle is rewritten to "The rule" before any arm is run, and the figures
  of the arms with a default at K = 3 or with rows of weight zero are not carried.
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
- Consumers. bartCause's `bcf()` (dbarts-1.0, 6c1bff9; read again at 70ad798, the sd unit's edit, where
  `sd.control` and `sd.moderate` default to `NULL`) passes `amplitude.prior.variance = b.prior.variance`
  on its treatment forest and holds either forest with `amplitude = fixed()`; it reads
  `getLeafPrior(1L)$response.scale` and `$response.shift` and the coefficients as `glue`, and nothing else of
  this surface. Two of its test files build the same sampler by hand with `amplitude.prior.variance = 0.5`.
  Its treatment basis is a two-level factor from the kind-by-class slice on. stan4bart (bartcore,
  a9d081b), treatSens (dbarts-1.0, aecec71) and bairrtt (main, 3f57f61) declare one forest with
  `dbartsForests$forest(n.trees = )` and read no coefficient, no forest's sd and no leaf prior of a model
  of several forests (searched again 2026-10-07).

## The rule

s is a forest's `sd`, in the response's units; U is the unit the sd unit slice records (sd(y - offset) over
the rows in the likelihood at creation under gaussian, 1 under probit, pi / sqrt(3) under logistic). Term f at row i is
(B_f(i, .) . a_f) F_f(x_i): the basis row times the forest's coefficients, times the forest's own fit. The
kind of a forest is what forest-kind-by-class records and fixes for the sampler's life (dec-B452):
nothing, the two levels of a factor, or a numeric column. Under dec-B354 a held coefficient is 1 and
the basis carries the coding; the engine keeps a two-level factor as two columns, on which that is the
pair (0, 1).

| the forest multiplies | coefficients | drawn | held, `amplitude = fixed()` |
|---|---|---|---|
| nothing | one | F has sd U; the size of a is half-Cauchy with median s / U: a F has median size s | a is 1; F has sd s, exactly |
| the two levels of a factor | one a level | each a_l is N(0, 1/2); F has sd s / 0.674: (a_2 - a_1) F has median size s | a is (0, 1); F has sd s, exactly, the second level against the first |
| a numeric column w | one | a is N(0, 1); F has sd s / 0.674: a F has median size s, per unit of w | a is 1; F has sd s, exactly, per unit of w |

1. Stated. `sd = s` is one number, per unit of the column as it was given where the forest multiplies a
   number, and nothing is read from the data.
2. Not stated. No basis, drawn: 2 U under gaussian, U under probit and logistic. No basis, held: under
   gaussian the size `bart()` gives the one forest of an ordinary fit on the same response, a quarter
   of the response's range as the engine holds it (about 1.4 sd(y) for a normal sample of 300), taken
   from the code that gives `bart()` its own so that the two cannot drift (dec-B290); U under probit
   and logistic (dec-B293). A factor: U. A number: s = U / sd(w). For a forest with a basis held or
   drawn alike, and whatever the number of forests (dec-B291).
3. The scale, sd(w): the sample standard deviation of the column, n - 1 in the divisor, unweighted, over
   the rows in the likelihood when the sampler is first created (dec-B303): the rows left by `subset`
   and the na.action that have a positive weight and that no mask leaves out. No later change of the
   weights or the mask moves it. Where creation read the response's scale from every row, no row
   being in or the rows in holding fewer than two distinct values of the response (dec-B302,
   dec-B309), the column's scale is read from every row too, under that one warning; under probit and
   logistic, where nothing is read from the response, that is where no row is in, and the warning is
   raised for the column ([Settled calls](#settled-calls), 1). The engine takes it, once, and
   R records it with the data; every later construction is handed the record and takes nothing again. A
   column constant over the rows read has no scale, and neither has one whose standard deviation is
   under 2^-26 of its largest absolute value there: where a default would need one it is refused by
   name, though the column may vary over the rows left out; with an sd stated it is created
   (dec-B307).
4. The column is never centred and never divided: the engine multiplies by it as given, at the fit, after a
   swap and at new rows. A default fit on w and on c w is the same fit, for any c that leaves both
   representable: a row leaves a forest's update when its multiplier is at most 2^-26 of the largest
   absolute multiplier that forest has in the sweep, and for no smaller reason.
5. A held value goes by the forest's kind and never by where the forest stands. `fixed()` takes no
   value but 1: another coding is another basis, a numeric column held at 1, and another size is `sd`
   (dec-B354).
6. The variance of a coefficient is not an argument: `amplitude.prior.variance` goes.
7. A swap moves the basis and nothing else: no scale, no sd in force, no coefficient, no tree, and
   never the kind (dec-B452). A numeric column of zeros is taken at a swap (dec-B346) and derives
   nothing, as any column.
8. A write, `$setLeafPrior(forests = )`, states s, as forest-sd-unit has it; the record of a scale stays.
9. Where a drawn coefficient starts is not part of the rule. A forest with no basis or a factor starts where
   it starts today, which is what keeps its draws; a numeric coefficient starts at 1, where today one
   that stands second starts at 0.

Before and after, per path.

| path | before | after |
|---|---|---|
| creation at every door: a `forest()` term at `bart()` and `dbarts()`; a `forests` list on a formula, on a matrix, in `dbartsSpec()`; `dbartsData(bases = )` | the table of Context | the rule; a numeric forest that states no sd has its scale taken and recorded |
| creation from a data object that carries a scale (`sampler$data`) | | used, for a forest that states no sd, whatever the weights then in force; checked there and nowhere else |
| a declaration that gives a forest of such an object another basis | | that forest's record is dropped and taken from the new column |
| `$setForestBasis` | the row norm is taken again: the forest's own sd follows the new column | nothing but the basis moves |
| `$setLeafPrior(forests = )` | the number goes to one of two channels | the forest is stated at s |
| `copy()`, a reload, `new("dbartsSampler", control, model, data)` | built from the column then in force | built from the record: the creation's prior, bit for bit |
| `setState`, a state in a reload | carries no prior | unchanged; a held block in the state is passed by, as today |
| `predict` at new rows, `fitted`, `extract` of draws | coefficients times the basis as given | unchanged |
| `sampleTreesFromPrior`, `sampleLeafParametersFromPrior` | each forest's own sd; coefficients left | unchanged, at the rule's own sd |
| `samplePriorPredictive` | stops inside `predict` | refused by name |
| `setData`; a warm start | refused | refused |
| `setResponse`, `setOffset`, `setWeights`, `setPredictor`, `setActiveRows` | move no forest's sd; with `updateScale = TRUE` refused on several forests (dec-B362, dec-A207) | unchanged, the refusal included; no change of weights or mask moves a recorded scale |
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
forest(x1, basis = factor(z), amplitude = fixed())` (K = 4). Under dec-B291 a default does not depend
on K, so a forest has the same figure in the table and in the model of four. Each literal a test pins
is computed in the test from its fixture.

| forest | `leaf.prior` | `sd` | `sd.stated` | `multiplier` | `amplitude` | `prior.sd.of` |
|---|---|---|---|---|---|---|
| no basis, default | `forest()` | `2.904` | `FALSE` | "none" | absent | "forest (prior median)" |
| no basis, held, default | `forest()` | a quarter of the response's range (dec-B290) | `FALSE` | "none" | `fixed()` | "forest" |
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
  forest's own sd, k being 1 on every forest of a model of several: [Settled calls](#settled-calls), 2),
  `response.scale`, `response.shift`. Gone: `amplitude.prior.variance` (push 1);
  `amplitude.prior.scale`, `leaf.scale.factor`, `leaf.scale.divisor`, `basis.row.norm` (push 2).
- `extract(fit, type = "leaf.prior.sd")` on that fit of four forests:
  `list(forest1 = 2.904, dose = c(dose = 71.18), age = c(age = 0.0827), "factor(z)" = 1.452)`; with
  `forest = "dose"` it is `c(dose = 71.18)`, the entry itself. `summary()`'s line under its table
  printed the same entries by label, as the sd unit slice built it, a numeric forest's number per unit
  of its column; after push 3 it does so only where no block prints, a fit made with `keepFits = FALSE`.
- The printed block (push 3, trimmed in its fix round: the package owner's ruling of 2026-10-22). `print`
  of a fit and `show` of a sampler give the short block, one line for each forest, with no range, no "per
  sd", no source and no mark of a held coefficient, under a plain header and over a last line that says
  where the full block is:

      forests:
        1 forest1: no basis; sd 2.904
        2 dose: numeric basis, mean 0.01994, sd 0.0204; sd 71.18 per unit
        3 age: numeric basis, mean 50.88, sd 17.55; sd 0.0827 per unit
        4 factor(z): factor basis, 2 levels; sd 1.452 between the two levels
      (each forest's prior in full: summary())

  A sampler's last line points to `$getLeafPrior()`, a sampler having no `summary`. `summary()` of a fit,
  after its table, and a sampler created with `verbose = TRUE` give the full block:

      forests, sd in the response's units (sd(y) 1.452):
        1 forest1: no basis; sd 2.904 (default)
        2 dose: numeric basis, mean 0.01994, sd 0.0204, range 0 to 0.0698; sd 71.18 per unit, 1.452 per sd (default)
        3 age: numeric basis, mean 50.88, sd 17.55, range 19 to 88.4; sd 0.0827 per unit, 1.452 per sd (default)
        4 factor(z): factor basis, 2 levels; sd 1.452 between the two levels (default); coefficient held

  A stated sd reads `sd 30 per unit, 0.612 per sd (stated)`. Under probit the header is `forests, sd on the
  latent scale (error sd 1):` and under logistic `(error sd 1.814)`. A column's mean, sd and range are
  those of the column in force over the rows in the likelihood when the block is printed; where a
  default was taken from another column, or from other rows, the full line ends `(default,
  taken from a column of sd 0.0204)`. A constant column, which has a stated sd, prints no "per sd".
  A mean is printed through `zapsmall` against the sd. This is the
  one sign a 0/1 number (dec-B263: at 20 percent treated beside one other forest, `mean 0.1967, sd
  0.3981, range 0 to 1; sd 3.647 per unit` where the factor's line says 1.452, at any number of
  forests), a column far from zero
  (dec-B264: `mean 50.88, sd 17.55`), one wild value (the range of the full block, the mean and sd of
  the short) and a stated sd in the wrong units get: nothing warns. For a default, "per sd" is the same
  number on every numeric line of a model; it is there for the stated case.

## Refused forms, with their texts

Base R's style, as in written-surface. `<f>` is `forest 2`, or `forest 2 ("dose")` where the forest has a
label. The push that adds each is in brackets.

    [1] 'amplitude = fixed(2)': a held coefficient is 1, and what it multiplies is the basis: nothing, or 0 and 1 for the two levels of a factor; fixed() takes no other value. Write fixed(), and state the forest's size with 'sd'
    [2] 'amplitude = fixed(2)': a held coefficient is 1, and what it multiplies is the basis: nothing, a numeric column as given, or 0 and 1 for the two levels of a factor; fixed() takes no other value. Write fixed(), state the forest's size with 'sd', and write another coding as the basis, as basis = z - 0.5
    [2] <f> states no sd and its basis is constant (61.7 on every row in the likelihood), so there is no standard deviation to take a default from; state one with sd =
    [2] <f> states no sd and its basis is constant to rounding (its standard deviation is 1e-11 of its largest value), so there is no standard deviation to take a default from; state one with sd =
    [2] no observation is in the likelihood (0 of 150 have a positive weight), so the scale of <f>'s basis is taken from all 150            (a warning; where the response's own warning is raised, one warning names both)
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
"'amplitude.prior.variance' is the prior on a basis forest's amplitudes". Gone with push 2: the
refusal of a held coefficient on one numeric column (dec-A171, confirmed as dec-B353 until this is
built). The [1] text replaces the one [`validateForestAmplitude`](../../R/model.R) has today, which
names the same values and does not say that the basis carries them; dec-B354 has the refusal point to
the basis for a coding and to `sd` for a size, and the pointer to the basis waits for push 2, a held
numeric column being refused until then. Kept as they are: every text of the kind-by-class slice, whose refusals
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
    A number w. A column that is an ordinary covariate should be centered: basis = scale(w), or
    basis = I(w - c) for a reference value c. The other forests describe the response where w is 0, and
    a column that never comes near 0, an age or a year, makes them an extrapolation and the fit poor,
    with nothing in the chain to show it. A dose or an exposure, where 0 means none, is left as it is.
    The forest enters as w (a f(x)), and sd is the size of a f(x), the change in the response for one
    unit of w; the term at an observation has size |w| sd. The column is used as given: it is never
    centered and never divided, so a stated sd is per unit of w, and the same sd with w in grams and in
    milligrams states two priors a thousand times apart.

    Where the coefficient is held (amplitude = fixed()) it is 1, what it multiplies is the basis, and
    "size" is the prior standard deviation, exactly. Where it is drawn the size is itself uncertain,
    and sd is its prior median.

    Not stated, sd for a forest with no basis is 2 sd(y) under a continuous response where its
    coefficient is drawn, the size bart() gives its forest (a quarter of the response's range) where
    it is held, and one unit of the link's error under "probit" and "logistic" either way; sd(y), or
    that unit, for a factor; and for a number that same size for one standard deviation of the column,
    sd(y) / sd(w) per unit, held or drawn. None of these depends on how many forests the model has. A
    default fit is then the same fit whatever units w is in. sd(y) and sd(w) are taken once, over the
    observations in the likelihood when the sampler is created (a positive weight, and not masked
    out), and kept: replacing the column, the weights or the mask later does not move the prior. A
    caller who replaces a column and wants the default of the new one states it,
        p <- sampler$getLeafPrior(2)
        sampler$setLeafPrior(forests = list(forest(), forest(sd = p$sd * p$basis.scale / sd(w.new))))
    A column that is constant over those observations has no standard deviation and needs sd stated.

    A treatment belongs in a factor or a logical vector, where sd is the treatment effect. Written as a
    0/1 number it is scaled as any number: sd(z) is 0.5, 0.4 and 0.22 at 50, 20 and 5 percent treated, so
    the default effect is 2, 2.5 and 4.6 times the one factor(z) gets;
    sd = sd(y) on the number states the factor's. The same holds for a coding written as a number and
    held, basis = z - 0.5 with amplitude = fixed().

    A held coefficient cannot adapt to an sd that is too small. At the default it is sized for a slope of
    about one sd(y) per sd(w); where the effect may be larger, under "probit" and "logistic" above all,
    state sd or let the coefficient be drawn.

After push 1 the item has the paragraphs on a held coefficient and on a forest with no basis and a
factor, with the defaults of dec-B290 and dec-B291, and says of a number what the sd unit slice's
help says with its default no longer divided by the number of forests. Until push 3 it does not say that `print`
shows what each forest resolved to; push 3 adds the sentence.

## Constraints

- A fit whose forests have no basis or a factor, coefficients drawn, draws what it draws on the base build,
  to the bit, in every push, stated sds included, at two forests; at three or more it does so from push
  1's first commit on, and across that commit where every forest with a basis states its sd, or
  against the base build stating the unit (the `rng:` line's identity). In push 1 a fit with a drawn
  numeric multiplier does too, on the same terms.
- The help promises a default fit that does not depend on a column's units only on a tip that has the
  relative bound: the bound is push 2's first commit and the sentence its last.
- The law is resolved in one function of the engine and nowhere in R: R states what the caller stated and
  the kind, and reports what the engine returns. The engine holds every default, the column's scale
  among them: R computes no standard deviation that the engine uses.
- A forest's law is resolved once for each construction, by the chain, and the combiner is handed the
  result; no second derivation can differ from the first.
- "Constant" is decided by comparing values, never by a computed standard deviation being zero, and in
  the engine, which knows the rows the scale is read from; R gives the engine's refusal its text.
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
- The kind is read from the one function forest-kind-by-class added,
  [`forestBasisKind`](../../R/forestBasis.R), and its swap table stands as it is:
  no forest changes width or kind (dec-B452), so this slice adds no row to it, and the kind stated to
  the engine at a construction holds for the sampler's life.
- The rows a column's scale is read from are asked of the engine once, where the unit's are; R holds
  no second test of which rows are in.
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

1. Defaults and held coefficients. First commit, gated by its identity: the default of a forest with a
   basis is one unit whatever the number of forests (dec-B291). Then the holds: the engine is told
   each forest's kind; the held values go by kind, a held forest's
   sd is exact, a held forest with no basis that states none takes `bart()`'s size under gaussian
   (dec-B290), forest-defaults-by-kind's two interim refusals are lifted, and
   `amplitude.prior.variance` goes. The unit on a constant response is the size one forest gets
   (TODO constant-response-copy-refused, closed by the push). The numeric rule is otherwise untouched, so dec-B353's refusal of a
   held numeric column stays. The exact gate's script is created with its held arms. bartCause's edit the
   same day.
2. The numeric multiplier. First commit, gated alone and bit for bit: the bound under which a row leaves
   a forest's update becomes relative. Then the law: the column as given, the normal coefficient, the
   default for one standard deviation of the column, over the rows in the likelihood, with its record,
   a held numeric column, the reader's and `extract`'s last change. dec-B353's refusal is lifted. The gate's numeric arms and four
   baseline scenarios. The design note is completed here.
3. The printed block. R only.

No tip between them fits a model nobody asked for: after push 1 a held numeric column is still refused,
not held under half a law; push 2 lands the numeric law, its reader and its help together. Push 2's tip
lacks the printed block for as long as push 3 takes; its help does not promise one. Push 3 changes no
file under src/, so it is built in its own worktree while push 2's engine runs are watched, and takes
the R-surface tier's review.

## Steps

"Fails today" is what the base build does where the test expects otherwise. New names are for the reader's
sake. Calls are in push 3 of written-surface's spelling; sds are in the response's units. Fixture, unless
said: 150 rows with x1, x2, a 0/1 z, its factor zf, dose (sd 0.02), age (mean 50, sd 17), an offset column
and weights of which ten are 0; a response whose standard deviation is neither 1 nor its range.

### Push 1: held coefficients

1.0 First commit: the default of a forest with a basis (dec-B291).
    [`defaultForestSd`](../../src/bartcore/combiner.hpp) gives 1 for a forest with a basis whatever
    the number of forests, and no longer takes the count. Nothing else of the engine moves in this
    commit. Tests: tests/cpp, the calibration map's literals at three and four forests, a factor and a
    number, each one unit (fails today: sqrt(2 / 3) and sqrt(2 / 4)), and at two forests every leaf
    scale the base build's to the bit; tinytest, in the sd unit slice's file, the reader's `sd` of a
    default factor forest at two, three and four forests equal to the unit, a default numeric
    forest's the same at each count, and under gaussian, probit and logistic a model of three forests
    `identical()` over 20 sweeps to its twin whose forests with a basis each state the number the
    reader gives. Across builds the `rng:` line's identity is a row of the pair script (Verification):
    the new default against the base build stating the unit. The
    pins of sqrt(2 / K) in test-bcf-family.R, 16 assertions by the sd unit's recheck, and that
    slice's own at three forests take 1. man/forest.Rd's sentence on the default loses its K. Gated
    with the push; the three compares and the four snapshot files are run after this commit alone,
    every scenario identical, there being no recorded model of three forests.
1.1 The engine. [`ForestSpec`](../../src/bartcore/combiner.hpp) states a kind (nothing, levels,
    numbers) beside the [`ForestStatement`](../../src/bartcore/combiner.hpp) the sd unit slice has it
    carry, in place of the one bit it has today:
    retired: [`ForestSpec::createdWithoutBasis`](../../src/bartcore/combiner.hpp), now gone, and the forests' record's
    `hasBasis` go from what the bridge reads, no forest gaining a basis after creation (dec-B452).
    [`expandForestSpecs`](../../src/bartcore/combiner.hpp) states BCF's two forests as nothing and
    levels. [`forestMapEntry`](../../src/bartcore/combiner.hpp), which holds the rule already, takes
    the kind where it takes that bit, and the held rows of "The rule" for a forest with no basis
    and a two-level block: a held forest with no basis has own sd s, by the leaf scale and not the
    half-Cauchy median, and where it states none under gaussian s is a quarter of the response's
    range, by the expression that sizes the one forest of a single-forest chain at k = 2, called and
    not copied (dec-B290; under probit and logistic one unit, dec-B293); a held two-level block has no
    divisor. A held block is set to its kind's value
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
    from the data object, as [`forestBasisKind`](../../R/forestBasis.R) reads it: no basis; a basis
    with the levels forest-kind-by-class records (`basis.levels`, R's
    [`dataBasisLevels`](../../R/data.R) answering `NULL` for an object saved without the slot); a
    basis without.
    [`applyForestStatements`](../../src/R_interface_bartcore.cpp) no longer reads a variance or
    `hasBasis`.
    [`bartcore_getLeafPrior`](../../src/R_interface_bartcore.cpp) reports the kind the engine holds and
    the hold. Tests: no loop of its own. forest-kind-by-class's loop over the creation doors
    (["The seven creation doors"](../../inst/tinytest/test-forest-kind.R)) and its test of
    ["the same sampler again"](../../inst/tinytest/test-forest-kind.R) (a copy, a reload,
    `new("dbartsSampler", control, model, data)`) each gain one assertion a cell: the kind the bridge
    reports is the kind [`forestBasisKind`](../../R/forestBasis.R) gives for the data object (a kind
    dropped between the data and the engine reads "numeric" for a factor at that door).
1.3 R. [`forest`](../../R/model.R) loses `amplitude.prior.variance`, and
    [`validateForestKnobs`](../../R/model.R), [`resolveForests`](../../R/model.R) and the multi-forest
    block of [`resolveSamplerSpec`](../../R/spec.R) their lines for it.
    retired: [`refuseHeldShape`](../../R/model.R), which forest-kind-by-class left with three rows (one numeric
    column; a forest with no basis second; two columns not second), keeps one,
    dec-B353's, read from the forest's kind and no longer from a width, and loses the two by position:
    a forest with no basis and a two-level factor are held wherever they stand.
    [`validateForestAmplitude`](../../R/model.R), which refuses a value other than 1 already, takes
    the [1] text. The reader
    ([`reportLeafPrior`](../../R/dbarts.R)) gains `multiplier` and `amplitude` (`fixed()` or absent),
    takes the four `prior.sd.of` values of a forest with no basis or a factor, and loses
    `amplitude.prior.variance`.
    Tests, test-forest-held.R:
    - By kind, at every door and position. `amplitude = fixed()` on a forest with no basis, on
      `factor(z)` and on a logical, in a formula term, a `forests` list on a formula and on a matrix,
      `dbartsSpec()` and a data object's `bases`; as the first, second and third forest, and on a forest
      with no basis standing second in a `forests` list (a formula puts it first wherever its term is
      written): `getForestAmplitudes()` is 1, or (0, 1), before and after 50 sweeps
      (fails today: refused at the positions forest-defaults-by-kind closed; before that slice (1, 1)
      third and first, and 0 for a plain forest second).
    - Exact. The reader's `k.scale` of a held forest equals its `sd`: stated 0.3 and default on a forest
      with no basis, stated 0.6 and default on a two-level factor, under each family, to 1e-12 (fails
      today: the response's scale, and s / 0.674). The default of a held forest with no basis under
      gaussian is the `k.scale` over k of a `bart()` fit of one forest on the same response, to 1e-12,
      and with ten rows of weight zero far out too (dec-B290); under probit and logistic it is the
      unit. A held factor's default is the unit at two, three and four forests.
    - What a fit records. `fit$fixed$amplitude` has a 1 for each held forest and no entry for a drawn
      one, at every position (dec-B376; built for the shapes held today).
    - The size, from the engine's prior draws. Always on, one a family, 1500 sweeps and 10 percent: the
      second level's term of a held factor and a held forest with no basis have the stated sd. `at_home`,
      5000 sweeps a row and 4 percent: the same at three sds.
    - Refused: `amplitude.prior.variance` is R's unused-argument error at both doors and through
      `do.call`; `fixed(2)` with its text.
    - Not moved: a model with a numeric column drawn beside a held factor keeps the base build's row norm
      (`k.scale` of the numeric forest to the bit), and dec-B353's refusal stands at creation.
1.4 The gates.
    - A new script, benchmarks/R/forest-prior-exact.R, which push 2 extends: one binary predictor, sigma
      held, the leaves integrated out in closed form; no quadrature in this push's arms. Arms: a forest
      with no basis held at sds of 0.3 and 2 of the unit and at its default, a quarter of the range
      the oracle takes from the response it drew; a two-level factor held at 0.6, at 0.2 and at
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
    forest-prior-exact.R and of bcf-latent-exact.R pass, with their gaps.
1.6 Respell. `amplitude.prior.variance` on the test lines that state it (31 lines of 6 files at the tip,
    recounted after the kind-by-class slice's removals): where it states 0.5 the argument is dropped and
    the fit is the same; where it states another value the test pinned the tip's law and goes with it,
    one pin of R's error taking its place. Five benchmark scripts drop the argument.
    forest-defaults-by-kind's pins of its two interim texts by position become step 1.3's creations. The
    reader's pins of `prior.sd.of` and `amplitude.prior.variance` for a forest with no basis or a factor
    take the new entries. Run the suite first and repair what it shows.
1.7 Help and records. man/forest.Rd: the usage, the `amplitude.prior.variance` item removed, the
    `amplitude` item (the held values, two in this push, wherever the forest stands; that a held
    coefficient makes `sd` exact; that a held forest with no basis has the size `bart()` gives its
    forest under a continuous response, and one unit of the link's error otherwise), the `sd` item's
    sentences for a held forest and for the defaults, the Details paragraph on the
    budget. [`dbartsSampler$getLeafPrior`](../../man/dbartsSampler-class.Rd) with its docstring. A new
    docs/design/forest-multiplier.md with its index row: the rule's held half, the changed sequences
    with their oracles, the re-recorded scenario, what the engine keeps general. docs/design/bcf.md and
    multiplier-combiner.md where they give the 0.674 to a held forest or name the variance. TODO:
    `forest-prior-args`.
1.8 bartCause, same day (its own commit on dbarts-1.0; its basis has been `factor(z, levels = c(0,
    1))` since 3ecf9fb, the kind-by-class slice's edit). In R/bcf.R: `b.prior.variance = 0.5` leaves the
    formals of `fitBCF` and of `bcf`; `amplitude.prior.variance = b.prior.variance` leaves the treatment
    forest's call; `b.prior.variance = b.prior.variance` leaves the call of `fitBCF`; and beside its
    refusal of `forests` and `bases`, `bcf` refuses the argument by name ("'b.prior.variance' is no
    longer an argument: the treatment coefficients' prior variance is fixed; state the effect's size with
    'sd.moderate'"), since its dots would otherwise take it in silence. man/bcf.Rd: the usage line and
    the item go; `update.a, update.b` reads "when FALSE the coefficient is held, at 1 for the prognostic
    forest and at 0 and 1 for control and treated, and sd.control or sd.moderate is then that forest's
    prior standard deviation, exactly". tests/testthat/test-14-bcf.R and test-03-responseFit.R: the two
    hand-built samplers drop `amplitude.prior.variance = 0.5`, the same fit. Where `update.a` is `FALSE`
    and `sd.control` is not given, `bcf` states `sd = 2 * sd(y)` on its first forest, over the rows it
    hands dbarts with a positive weight, bcf's own convention, and its help says so (dec-B290). A
    `bcf()` fit with
    `update.a = FALSE` or `update.b = FALSE` moves; no test there pins a draw under either. stan4bart,
    treatSens and bairrtt are not edited, at this push or the next two.
1.9 The unit on a constant response (TODO constant-response-copy-refused, closed in this push). Under
    gaussian a response whose rows read hold one value, every row by dec-B309's fallback, has an
    anchor s of 0 at [`Chain::latentScaleAnchor`](../../src/bartcore/chain.hpp): every forest's
    default is then 0, a stated sd is refused for want of a unit, and `copy()` of the sampler is
    refused by the bridge's check of the recorded anchor. The coordinator's call (2026-10-17): the
    anchor on such a response is a quarter of the window of width one the response model gives it,
    the size the one forest of a single-forest chain has at its defaults, by the expression step 1.1
    calls for a held forest with no basis; the unit, s times the window, is then a quarter of it
    too, so forests of several take the sizes their defaults give against what one forest gets. Under
    probit and logistic the anchor is the link's and nothing moves. Nothing else of the anchor
    moves: a response with any spread over those rows keeps its bits. Tests, tests/cpp: a K-forest
    chain on a constant response has anchor and unit a quarter of the window, its leaf scales those
    the map gives on that unit, a stated sd is taken, and a response with spread keeps its anchor
    to the bit; tinytest, in test-forest-sd-unit.R, the TODO's reproduction: created with `sd` 0.5
    and 0.25 of the response's units for the forest with no basis and the factor, a stated sd taken,
    `copy()`, a reload and `new("dbartsSampler", control, model, data)` created with the creation's
    `sd` and `k.scale`, and the copy's next draws `identical()` to the original's; one forest on the
    same response unchanged (`k.scale` 0.5).
1.10 Mutations (Verification).

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
    over the rows it is handed, n - 1, the mean taken first and the deviations from it squared, as R's
    `sd()` does. The chain hands it the rows of rule 3 (dec-B303): the rows the response's range was
    read from ([`GaussianResponse::scaleRows`](../../src/bartcore/model.hpp), every row where it
    answers none), and under probit and logistic, which record none, the rows in the likelihood by the
    membership the response model already tests (a positive weight and no mask), every row where none
    is in; the chain reports that fallback as a status, the engine never warning, and the R interface
    raises the one warning. Nothing reads the weights' sizes.
    It throws for a column with no scale: fewer than two rows; every value equal to the first; a
    standard deviation at most 2^-26 of the largest absolute value. Tests, tests/cpp: columns of sds
    0.25 and 8 give exactly those; a 0/1 column at 20 percent; a column with half its rows zero equals
    the literal over all the rows handed; with a fifth of the rows at weight zero and values far out
    on them, the literal over the rows of positive weight (fails over every row), under gaussian and
    under probit; with every weight zero, the literal over every row and the status; with one row in
    under gaussian, the same (the response's own fallback), and under probit the throw for a constant
    column; a column that varies over every row and is constant over the rows in throws; 1e7 + i for
    i in 1 to 150 equals the literal R gives to 1e-10 relative
    (fails with one pass); the constants 1, 0.5, 0.1, 61.7 and 1e8 + 0.1 each throw (the last three fail
    with a test of the computed sd against zero); 61.7 + 1e-9 sin(i) throws; one row throws.
2.3 The law. [`ForestSpec`](../../src/bartcore/combiner.hpp) gains a scale (`basisScale`, not a number
    for "take it"). The numeric row of [`forestMapEntry`](../../src/bartcore/combiner.hpp), which takes
    the scale where it takes the row norm today, is "The rule"'s: s stated, or the default over the scale
    handed or taken; the forest's own sd s / 0.674 where the coefficient is drawn, under variance 1 and
    starting at 1, and s where it is held at 1. A numeric block of several columns throws, and a
    fixture of tests/cpp that builds one is restated on one column, or as a levels block where that is
    what it stood for.
    retired: [`basisRowNorm`](../../src/bartcore/chain.hpp) goes with everything that read it. One function of
    the chain sets a forest's leaf scale from the sd in force,
    [`mapLeafScale`](../../src/bartcore/chain.hpp) as the sd unit slice has it for a write and a swap,
    and the constructor calls it too, where today it writes the expression out.
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
      creation: the scale is `sd(dose[kept])` to 1e-12, `kept` written in the test, leaving out
      the rows of weight zero, whose dose is far out, and ignoring the later mask; and a default
      forest's `sd` is `U / sd(dose[kept])`, U the sd of the response over the same rows, at two and
      at three forests. After `setWeights` brings the rows of weight zero in, and after
      `setActiveRows`, the record and `sd` are `identical()` to creation's. Created with every
      weight zero: one warning, and the scale is `sd(dose)` over every row. A dose constant over the
      rows of positive weight that varies over the others: refused with no sd stated, created with
      one. Under probit the same four, the warning included.
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
2.6 Refusals and what is accepted. The engine refuses, at creation, a numeric column that is constant
    over the rows its scale is read from, by comparing its values, or whose standard deviation there is
    at most 2^-26 of its largest absolute value, on a forest that states no sd; the R interface gives
    the refusal its [2] text, with the forest's label and the constant's value, from what the engine
    reports, so that which rows are in is decided in one place. With an sd stated the column is
    created (dec-B307). [`samplePriorPredictive`](../../R/dbarts.R) refuses several
    forests at its door. dec-B353's refusal goes, and
    retired: [`refuseHeldOneColumn`](../../R/model.R) and what is left of
    retired: [`refuseHeldShape`](../../R/model.R) with it: a numeric column is held at 1.
    Tests: `rep(1, n)`, 0.1, 61.7 and 1e8 + 0.1 on every row, and a column constant on the rows `subset`
    keeps, each refused at four doors, in the [2] text with the forest's label, and created with
    `sd = 2`; the same five through
    `new("dbartsSampler", control, model, data)` with the record cleared (fails on a check of the
    computed sd: the last three are created);
    61.7 + 1e-9 sin(i) refused, and 1e7 + i created with the scale R's `sd()` gives; the tests of
    dec-B353's refusal turned over: `amplitude = fixed()` on dose is created at every door, held at 1 as
    the second, third and first forest, with `k.scale` equal to the sd (stated, and
    `U / sd(dose)`) and a 1 in `fit$fixed$amplitude`; swapped to another column, a column of zeros
    among them (dec-B346), it is still held at 1 with the same
    `k.scale`; a coding held as a number, `basis = z - 0.5`, is created and held at 1, and its term
    is the forest's fit times -0.5 and 0.5, by `type = "forest"`, to 1e-12; a held 0/1 number and the held factor of the same column, at one stated sd, have
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
    - forest-prior-exact.R gains the numeric arms, nine of one column (Context) and the arm of three trees
      with a numeric column, drawn and held, in the place of its second forest; the oracle written from
      "The rule" with R's `sd()` over the rows of positive weight. Functionals: the posterior mean of
      the term in each cell at each value of the basis; the probability that the multiplied tree
      splits; the probability that the coefficient is within 0.5, 1 and 2. The verdict and the seeds are
      push 1's. What it does not see is said in its header: n for n - 1 in the scale on a held arm
      (step 2.2's literals and step 2.5's 1e-12 hold that), a numeric multiplier under logistic, and
      `subset`.
    - The prior from the engine, tinytest. Always on, one a family, 1500 sweeps and 10 percent: a
      default numeric forest has size U per standard deviation of its column, in a model of three
      forests. `at_home`,
      5000 sweeps a row and 6 percent: default and stated, a column times 1000, an age, a 0/1 number at
      20 percent, held, under each family.
    - The units twin, tinytest, always on: dose against 1024 dose and against 2^-40 dose, default, drawn
      and held, gaussian and probit: identical draws.
2.9 The baseline. Four scenarios are added to [bcf-equivalence.R](../../benchmarks/R/bcf-equivalence.R)
    and recorded: a numeric column at its default; a numeric column through a swap, a write and a
    copy; a held numeric column; a held forest with no basis under gaussian. The 15 in force are
    bitwise on the push's build first, counted; the file then holds 19. The MANIFEST row names the
    oracle of the four: forest-prior-exact.R's arms with their gaps, and the units twin. They are there
    for what comes after: the slices that move the forests' record and fold `basis.scale` into the
    data's other records are gated bit for bit, and until now no recorded draw would move if one of
    them took a scale twice or dropped a record.
2.10 Respell and repair. Run the suite first. Expected: the assertions that pin the row norm and a
    numeric coefficient's variance (retired: ["medianRowNorm"](../../inst/tinytest/test-bcf-family.R), 13 at the
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
    its column's name. docs/design/forest-multiplier.md completed: the rule, the scale and its
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
  exit status; in pushes 1 and 2 the four BCF scripts and forest-prior-exact.R also in full.
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
- A hot path is touched once, by push 2's first commit, which was compared against its parent:
  `bench-sampler.R record` on 738b1015 then `compare` on a51ac7fe, run by the orchestrator on the idle x86
  bench box (4 cores, load about 1), one machine: no metric regressed more than 5%. The ratios: embedded
  offset gibbs step 1.000, binary run 0.988, run t200 1.000, run t75 1.000, run n10000 0.994,
  setPredictor accept 1.035, setPredictor reject 1.000.

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
| 29 | three and four forests, a factor and a number at their defaults, drawn, gaussian, probit and logistic, against the base build with each of those forests stating the unit | identical, from the first commit | identical | the factor's forests identical; a model with the number differs |

Rows 03, 15 and 17 have three or four forests: where a forest with a basis is left at its default they
differ across push 1's first commit, and are identical in the form of row 29, the base build stating
the unit; with every such sd stated they are identical as the table says. A held forest with no basis
at its gaussian default (dec-B290) has no spelling on the base build, which ignores a held forest's sd:
its oracles are the exact gate's arm and step 1.3's equality with `bart()`'s own forest.
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
  forest-prior-exact.R (0.10 and 0.21 at its two stated arms on the draft);
- push 1, a held forest with no basis left at the response's scale: the same three (0.033);
- push 1, a held forest with no basis given one unit by default under gaussian: step 1.3's default;
- push 1, a leaf's sd not divided by the root of the tree count: the gate's arm of three trees (167
  on the draft);
- push 1, sqrt(2 / K) left in the default of a forest with a basis: step 1.0's literals at three and
  four forests, and the gate's arm of three trees, K = 3;
- push 1, a held forest with no basis left at 2 sd(y) by default, or given a quarter of the range of
  every row where rows are out: step 1.3's equality with `bart()`'s forest, and the gate's default arm;
- push 1, `amplitude.prior.variance` taken and ignored: step 1.3's unused-argument pin;
- push 1, a write to a held forest with no basis sent to the half-Cauchy median: step 1.1's write;
- push 1, a held block installed from a state: step 1.1's state;
- push 1, one of forest-defaults-by-kind's interim refusals left standing: step 1.3's creations;
- push 1, the anchor of a constant response left at 0: step 1.9's tests/cpp anchor and tinytest copy;
- push 2, the bound left absolute: step 2.1's twin, and step 2.8's units twin downward;
- push 2, the bound taken relative to the largest multiplier of every forest, or of the last sweep: step
  2.1's row at 2^-30 of its own forest's largest;
- push 2, a zero multiplier kept in the update: step 2.1's held block;
- push 2, the default per unit (the scale read as 1): step 2.5's `U / sd(dose[kept])`, and
  the gate (0.06 to 0.85);
- push 2, a scale 3 percent off: the gate, every arm (4.6 to 10.6 on the draft);
- push 2, n for n - 1: step 2.2's 0.25 and 8, and step 2.5's 1e-12 (the gate sees it on the drawn arms
  only just);
- push 2, the scale over the nonzero rows; over every row where rows are out of the likelihood; before
  `subset`; with weights; over the rows a later mask leaves; the fallback to every row dropped, or
  taken where the rows in vary: step 2.2's half-zero column and its rows of weight zero, step 2.5's
  "which rows", and the gate's weights arms;
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
| | a record a call can rewrite, and the arithmetic to restate a forest from a scale | the slot of step 2.5; `basisColumnScale` of step 2.2, which takes the rows it is handed and so the rows in at a later call (dec-B303), and the chain's one expression for a leaf scale of step 2.3 |
| | no state, no shipped signature: the update is a virtual of the engine's own facade and one more bridge argument | Constraints |
| `leaf.prior = normal(sd = )` on `forest()`, and `normal()` for the default back (dec-B266, dec-B276) | `leaf.prior` on `forest()` is R's unused-argument error | written-surface's pin |
| | the reader already says "stated" or "the default" | forest-sd-unit's `sd.stated`; no test here pins the class of `leaf.prior` |
| | a stated forest keeps its record, so a default put back has the creation's scale | step 2.5's "stated, then put back" |
| a heavy-tailed coefficient, `amplitude = student(3)` (dec-B262, dec-B271) | `amplitude` takes `fixed()` and nothing else | written-surface's pins |
| | one row of `forestMapEntry`: a drawn numeric or factor block under the forest with no basis's law | the law is one function |
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

- On push 3 of written-surface, landed: a basis's column name, which names a numeric forest's `sd`; the labels
  that name the reader's and `extract`'s lists and the texts; every call here is in its spelling.
- On forest-defaults-by-kind, landed: a forest with no basis may stand anywhere in a `forests` list, so
  step 1.3 holds one standing second; its interim refusals of held shapes, two of which push 1 lifts;
  selection by label in the reader; its line of tree counts, which the block follows.
- On the sd unit slice, built and landing first: the unit and the rows it is read from, the statement
  the engine is handed and [`forestMapEntry`](../../src/bartcore/combiner.hpp), which holds the
  rule, the four-number record, the reader's `leaf.prior`, `sd` and `sd.stated`, `extract`'s list by label
  and its entry for one forest, `summary()`'s line, a write that states. It left three things here by
  name: the default of dec-B291, which its open call gave to this slice; the default of dec-B290; and
  what k and `k.scale` are on a forest of several ([Settled calls](#settled-calls), 2). This slice changes
  rows of that function, puts the kind in the statement where one bit stands, and names one entry.
- On forest-kind-by-class, push 1 landed and push 2 built: the record of a factor's levels
  (`basis.levels` on the data object, read through [`dataBasisLevels`](../../R/data.R), whose
  writers and validity `basis.scale` is written beside) and the one function that gives a kind,
  [`forestBasisKind`](../../R/forestBasis.R), which the bridge's reading mirrors; its loop over the creation doors
  and its test of the same sampler again, which step 1.2 extends; a kind that no swap changes
  (dec-B452), which is what lets the engine be told once; its
  second push, after which every forest is one numeric column or two levels and none changes width, so
  that what moves in push 2 here is a number that was meant as one, the tests' blocks of indicator
  columns are the factors they stand for, a column of zeros is refused at creation at every door, and
  retired: [`refuseHeldShape`](../../R/model.R) has lost its row for three or more columns and
  retired: [`refuseHeldWidthChange`](../../R/model.R) is gone; bartCause's factor line. If that slice's build
  moves a name cited here, the names are read again before push 1 is built.
- To recheck once each has landed: the stand-in is of no further use for counts, most of what it failed
  being repaired or removed by then; run the suite on the push's first build and count. The pair
  script's rows. That no baseline scenario has gained a number or a hold. That `samplePriorPredictive`
  still cannot predict a sampler of several forests (if it can, it returns draws with the coefficients
  at their starting values, and step 2.6's refusal moves to push 1). That the suite still holds only the
  two held shapes forest-defaults-by-kind kept. The names of the functions cited here.
- [leaf-conversions.md](leaf-conversions.md), cross-family-state-install and the removal of the
  state's family record (dec-B283) have landed, as have response-scale-rows through its slice C and
  the k-internal slice (dec-B416), which left several forests out of its scope: nothing of this slice
  waits on them, and none edits a line this slice does.
- With slice D of response-scale-rows (TODO, not planned; dec-B304): a linear or gp leaf's covariate
  center and scale over the rows in. It needs the engine's answer to which rows are in for a family
  that reads nothing from the response, which step 2.2 adds; whichever is built second uses the
  first's.
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
- `updateScale = TRUE` on a sampler of several forests: refused for 1.0-0 (dec-B362), on the four
  setters that take it (dec-A207).
- Several forests with no basis sharing one total (dec-B292, TODO `forests-shared-total`): after the
  merge. Until then a constant column with a stated sd is the one spelling of a second such forest
  (dec-B307), created here drawn or held.
- To TODO as new entries: a drawn coefficient starts by its forest's position (1, and 0 for the second
  forest's first), a start by kind being a change of draws for a model of three forests;
  `samplePriorPredictive` for several forests, which needs the coefficients drawn from their prior; a
  consumer that creates a sampler of several forests through the shipped header gets a scale taken at
  each creation, there being no R object to record it on (none does).

## Ruled after this plan was amended

Worked into the text above by the recheck of 2026-10-11; this section is kept as the list of what was
ruled between the amendment and the recheck, with where each now stands.

- dec-B303, with dec-B302 and dec-B309: rule 3, the help's `sd` item, steps 2.2, 2.5, 2.6 and 2.8, and
  settled call 1 for the one case the three leave in tension.

      forest(x, basis = dose)          # default sd per sd(dose[weights > 0]), held

- dec-B290 and dec-B293: rule 2, the reader's table, steps 1.1, 1.3, 1.4, 1.7 and 1.8.

      forest(x1 + x2, amplitude = fixed())            # own sd range(y) / 4, as bart()'s forest
      forest(x1 + x2, amplitude = fixed(), sd = 2 * sd(y))   # what bartCause's bcf passes

- dec-B291: rule 2, step 1.0, and every place that had sqrt(2 / K). The sd unit slice's open call left
  it to this slice: its tip still divides.
- dec-B353: the refusal of a held coefficient on one numeric column stands until push 2 and is cited
  by that number.
- dec-B354: rule 5, the two texts for `fixed(2)`, step 1.1's holds by kind.
- dec-B307, dec-B346, dec-B452: rules 3 and 7, the constraint on the kind, step 2.6.
- dec-B362, dec-A207, dec-B361, dec-B376, dec-B418: the path table, Context, steps 1.3 and 2.6.

## Settled calls

The recheck left two calls open. The coordinator settled both on 2026-10-17, each as its option A,
to which the steps were written, and told the maintainer; each goes to the ledger with the landing.

1. Which rows a numeric multiplier's default reads when too few rows are in the likelihood at
   creation (dec-B296, dec-B302, dec-B303, dec-B309). The column's scale is read over the rows the
   response's scale was read over: wherever creation took the response's scale from every row,
   no row being in or the rows in holding fewer than two distinct values of the response, the
   column's is taken from every row too, under the same one warning, which names the column beside
   the response; where the rows in gave the response a scale, the column is read over them and is
   refused when constant there, as dec-B303 has it. Under probit and logistic, where nothing is read
   from the response, the column reads the rows in, and every row with the warning where none is
   in; there a sampler created with exactly one row in and a default numeric forest is refused
   where the gaussian one is created with a warning. The reason is dec-B303's own: a prior that sizes
   the column times the forest against the response's spread is coherent only when both are
   measured over the same rows. Run on the sd unit slice's build: a gaussian sampler with one row in
   takes the response's scale from all 150 rows under the warning and is created, as with none in,
   while two rows in read the two. Rejected: dec-B303 to the letter, under which a gaussian sampler
   with one row in is created or refused by whether a forest multiplies a number at its default;
   and a fallback of the column's own, which reverses dec-B303's sentence on a column constant
   over the rows in. Built in push 2 (rule 3, steps 2.2, 2.5 and 2.6).
2. What k and `k.scale` are on a forest of a model of several forests (dec-B201, dec-B416). As
   built: k is 1 on every forest of such a model, `k.scale` that forest's own prior sd before its
   coefficient, and `extract(type = "k")` and `fit$fixed$k` 1 for each forest. What k means on
   several forests is decided together with forests-shared-total (dec-B292), before 1.0-0, since
   the forests that share one total are the ones for which k is a parameter again and must define
   it anyway; the entries are final before the release. Until then a held forest with no basis at
   its gaussian default, which has `bart()`'s size (dec-B290), reports k 1 against a `k.scale` of a
   quarter of the response's range, where `bart()` reports k 2 against half of it. Rejected:
   `k.scale` as the data's scale on every forest, which moves every fit of several forests in its
   last bits; and no k on several forests, against dec-B193.

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
- Of the critique's questions: (a) the default of a held forest with no basis was put and is ruled
  (dec-B290, dec-B293); (b) which shapes may be held is settled by dec-B281 and dec-B282;
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
- The scale is unweighted and is read over the rows the response's scale is read over (dec-B303 in
  place of this plan's "counts rows of weight zero"): one estimator for both, and the one R's `sd()`
  gives for the column over those rows.
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

The recheck's own, 2026-10-11. The first eight are internal; the last five change what a caller sees
in a small way that no ruling names, and go to the maintainer with the landing's minor calls.

- dec-B291 is built as push 1's first commit and not as a push: it is one line of the engine, push 1
  changes a default already, and a battery of its own would repeat push 1's. It is gated by an
  identity the base build can state, which the rest of the push has for two of its sequences only.
- The kind stated to the engine replaces retired: [`ForestSpec::createdWithoutBasis`](../../src/bartcore/combiner.hpp), removed by push 1,
  and the bridge's reading of `hasBasis`, in push 1: forest-kind-by-class left the choice here, and
  after dec-B452 the bit records a case that cannot arise. The forests' record keeps `hasBasis` where
  R reads it.
- Step 1.2 writes no loop over the creation doors: forest-kind-by-class's has every door and object,
  and one assertion a cell there cannot drift from it.
- Which rows are in is decided in the engine alone: it takes the scale, it refuses a constant column,
  and R words the refusal from what it reports. A second test in R of which rows are in would be the
  kind of copy the sd unit slice's constraints forbid for the unit.
- The fallback to every row is a status the chain reports and the R interface warns on, the engine
  never warning (the shape dec-B405 gave a kept response scale).
- [`forestMapEntry`](../../src/bartcore/combiner.hpp) and [`mapLeafScale`](../../src/bartcore/chain.hpp)
  keep the names the sd unit slice gave them; the plan's own names for them are dropped.
- The [1] text's pointer to the basis for another coding (dec-B354) waits for push 2: until then a
  held numeric column is refused, and a message that sends the caller to one would send them to a
  refusal.
- Three pushes stand (Recheck against bartcore, "Size").
- The printed block summarizes a column over the rows in the likelihood when it is printed, prints no
  "per sd" for a constant column, and keeps the unit's number in its header under probit and logistic.
  dec-B347 keeps a scale the family fixes out of `summary()`'s line of fixed parameters; here it is
  the unit of every number in the block, and under logistic 1.814 is not a number a caller knows.
- `summary()`'s line keeps the form the sd unit slice gave it; a numeric forest's number there is per
  unit of its column, and `print`'s block is where the units are said.
- The reader's new entries are named `multiplier` ("none", "factor", "numeric"), `amplitude` and
  `basis.scale`, with six values of `prior.sd.of`, and it loses the five that decomposed the tip's
  calibration (dec-A108, which the triage marked covered by dec-B141: the reader reports the prior
  alone). A factor forest's `sd` is one unnamed number, dec-B253's difference between two levels.
- `samplePriorPredictive` on a sampler of several forests is refused by name, where today it stops
  inside `predict`; drawing the coefficients from their prior, which dec-B361's rule would then
  cover, is a TODO entry and no part of this slice.
- A coding written as a numeric column and held takes a number's default, per standard deviation of
  the column (dec-B263 as dec-B354 extends it), and the help says how to state the factor's.

## Recheck against bartcore

On 2026-10-11, against 37b59ea0, the tip of the sd unit slice as landed (read at its re-gate tip on
6f5b9610, the same code with that slice's 21 commits). Read: the ledger through dec-A207 and dec-B453 for every entry on a forest's
coefficient, a hold, `fixed()`, a multiplier, a forest's sd and its unit, a basis's scale, the number
of forests, the rows a scale reads and what is printed; [forest-sd-unit.md](forest-sd-unit.md) with
its code in the engine, the bridge and R; [forest-kind-by-class.md](forest-kind-by-class.md) as
rewritten with dec-B452 and dec-B453; the root TODO. Run, on that tip's shipped build, four scripts of
small single-chain fits, and one on 0.9-34, which has no `forest()`, no basis and no reader. Nothing
was built, and no suite, gate or compare was run.

Rulings.

| ruling | where the plan stood | now |
|---|---|---|
| dec-B253, dec-B260 to dec-B264, dec-B266, dec-B269 with dec-B280 and dec-B282, dec-B271, dec-B272, dec-B275, dec-B276, dec-B281 | carried | carried |
| dec-B265 | carried; `updateBasisScale` after the merge | carried; dec-B303 adds that it reads the rows in at its call |
| dec-B354 | in the header only; the two texts for `fixed(2)` sent the caller to `sd` alone | rule 5 and the texts, which point to the basis for a coding; the refusal of a value other than 1 is built already, in other words |
| dec-B353 | cited as dec-A171's refusal | cited by its number; built, and lifted by push 2 |
| dec-B290, dec-B293 | appended; the steps had 2 sd(y) for a held forest with no basis | worked in |
| dec-B291 | appended; sqrt(2 / K) in the goal, the rule, the help and the steps | worked in as step 1.0; the sd unit slice left it here |
| dec-B303, dec-B302, dec-B309 | appended; rule 3 counted rows of weight zero | worked in, with call 1, settled since |
| dec-B307 | not cited; the plan agreed | cited |
| dec-B346, dec-B452 | not cited | rule 7 and the constraint on the kind |
| dec-B362, dec-A207 | the plan agreed for the two response setters | the four setters |
| dec-B361 | not cited | the refusal of `samplePriorPredictive` stands beside it |
| dec-B376 | not cited | a held coefficient in `fit$fixed`, pinned in steps 1.3 and 2.6 |
| dec-B418, dec-B283 | the state "records its family" | the record is gone; an install converts nothing |
| dec-B416, dec-B201 | not cited | call 2, settled since |
| dec-B347 | not cited | weighed against the block's header (the recheck's calls) |
| dec-B289, dec-B308, dec-B453 | forest-kind-by-class's | nothing here |

Built by the sd unit slice, and so not this plan's to build: the statement handed to the engine and
the one function that resolves it, the unit and its record, the reader's `leaf.prior`, `sd` and
`sd.stated`, `extract`'s list by label, `summary()`'s line, a named single sd, a write that states.
Left by it to this plan: the default of dec-B291, that of dec-B290, the kind in place of one bit, the
channels' entries in the reader, and k on several forests.

What was false of the tip, and is changed above: the header's "none of which has landed"; that R and
the combiner were the measured tip's; "the scale can be taken on any of nine paths", which are nine
routes from R to two places of the engine; a state's family record; leaf-conversions as work still to
come; `forestMapEntry` and the statement as names to be invented; the suite's counts, which are left
as history with the tip's file count beside them; the reader's entries, which the sd unit slice
extended, and `summary()`'s line, which it added; bartCause's tip. What was checked and stands: the
table of the rule today, row for row; a held value by position, in three places; the row norm at a
swap, a copy and for a constant column; `samplePriorPredictive`; the refusals of `updateScale`.

Shared with forest-kind-by-class. It builds the record and the function this plan reads a kind
from, the loop step 1.2 extends, the slot pattern `basis.scale` follows, and it removes
retired: [`refuseHeldWidthChange`](../../R/model.R), now gone, and one row of retired: [`refuseHeldShape`](../../R/model.R); it
leaves this plan the two rows by position, dec-B353's row, the engine's bit, the creations on a 0/1
number and on a constant column, and the reader's and the block's words for a kind. No step is
planned twice. What this plan needs of it: one size a forest, a kind no swap changes, and every block
that stands for a factor written as one.

Size. Planned lines go from ~5400 to ~5800, the forecast from an upper figure of 8900 to a range of
9000 to 12500, by the ratio the last five slices ran at (the `budget:` line). Three pushes are still
the cut. Push 1 is the larger for two rulings and is ordered so that its first commit, the default,
is proved by an identity before the holds are built on it. Folding the printed block into push 2 would
put an R-only change of the R-surface tier under the engine tier's review and gates, and it can be
built beside push 2 as it is. Splitting the default out as a push would add a battery for one line.

## Push 1 as built

Built 2026-10-17 on origin/bartcore at 92ae2714 (forest-kind-by-class landed) and rebased onto 33f9c299
after its fix round, not landed. Commits, in
order: the plan's settled calls and cites (8c04ac4e); the default of dec-B291, the first commit
(967085ad); the engine and the bridge, steps 1.1 and 1.2 (bc93a28d); R, the help of `forest` and of
`getLeafPrior`, test-forest-held.R and the kind tests' assertions, steps 1.3, 1.2's tests and part of
1.7 (d97f0e05); the respelled tests, 1.6 (4f80c85b); the unit on a constant response, 1.9
(47838362, the last commit that changes package code); the gates, 1.4 (047d7d26); the BCF baseline,
1.5 (c16a0b83); the design note and records, 1.7 (03808f4a and the commits after it). bartCause's
edit, 1.8, is 30cd6fc and the fix round's 0057d2a and 4a70e7d on a branch wt/forest-mult-p1 off dbarts-1.0 at 3ecf9fb,
not on dbarts-1.0.

Against the steps:

- 1.0 as planned: `defaultForestSd` and `forestMapEntry` lose the count; tests/cpp pins one unit at
  K = 2, 3 and 4 for a factor and a number and the K = 2 leaf scale to the bit; the sd unit tests pin
  one unit at 2 to 4 forests and a three-forest default `identical()` to its stated twin under each
  family; test-bcf-family.R's pins take 1, the index prior now pinned at 1.04912 sqrt(K) s. On a
  reference build of this commit alone: 55, 15 and 11 identical, the four snapshot files passing.
- 1.1 as planned, with the calls below. `ForestKind` (none, levels, numeric) replaces
  `createdWithoutBasis`; the held rows are in `forestMapEntry`; the combiner holds a block by kind in
  `holdAtKindValue`. `testHeldCoefficients` and the restated fixtures, the drawn numeric forest
  keeping the base build's expression to the bit.
- 1.2 as planned: `parseData` reads `basis.levels` behind `R_has_slot`; the reader's raw matrix gains
  `multiplier` and `amplitude.held`. One assertion a cell in the door loop (by `shapeOf`, so every door
  agrees with the first) and in the same-sampler-again test.
- 1.3 as planned, test-forest-held.R with 362 results.
- 1.4: forest-prior-exact.R with the eight arms of the step, 51 functionals, in quick 32 seeds of
  25000 in about 7 s. The four BCF gates and sbc.R restated; bcf-latent-exact.R run first with its
  oracle as it was, six arms failing (worst |z| 10.70, 9.59, 6.90, 5.67, 5.68, 5.34; three forests
  under probit 3.97 and mode 2b 1.72 and 2.66 inside), then restated, worst 2.66.
- 1.5 as planned; MANIFEST row of bcf-equivalence-47838362.
- 1.6 as planned; the stated-variance pins went with their argument, R's unused-argument error in
  test-forest-held.R taking their place, and forest-defaults-by-kind's interim refusals became
  holds.
- 1.7: the help with the code (d97f0e05), the design note and the other notes after the baseline.
- 1.8 as planned, with the call below on a binary response.
- 1.9 as the settled step: `Chain::latentScaleAnchor` returns `singleForestDefaultSize` where the sd
  it reads is not positive.

Calls made building, where the plan was silent or wrong:

- The forests' record drops `hasBasis` with the bridge's reading of it. The plan kept it "where R
  reads it", and nothing in R does; a field nothing reads would let a test flip it and see nothing.
- A held forest with no basis carries no half-Cauchy scale: it has no mixture to draw, the reader
  gives it no `amplitude.prior.scale`, and a state's variance is not installed into it (restoreGlue
  installs a variance only where a mixture is). No draw depends on it.
- The held default of a forest with no basis is `singleForestDefaultSize()`, the engine's own
  `SamplerOptions` node scale over its k, handed to `forestMapEntry` in units of the anchor; the same
  function gives the anchor on a constant response, so the two cannot drift.
- Until push 2 removes `amplitude.prior.scale`, `$amplitude` on a drawn forest with no basis matches
  `amplitude.prior.scale` by prefix; the tests read `[["amplitude"]]`.
- The engine's refusal of a stated sd with no unit (`ForestSdUnitError` and the bridge's text) is no
  longer reachable from R, a constant response having a unit; it stays as the guard of a handed unit
  of 0, which tests/cpp pins.
- tests/cpp fixtures that stated a coefficient variance (1, 0.75) hand it to the combiner beside the
  map's priors; the K = 4 calibration fixture's variances are the map's 1/2.
- bartCause states 2 sd(y) on a held prognostic forest only where the family is not probit or
  logistic, read from the family handed over and otherwise from bcf's rule that a 0/1 response is
  probit (the fix round; it had read the response's values); under a binary family dbarts's default
  of one unit stands (dec-B293). Its test file gains a test of the refusal and of the held default.
- The "Out of scope" entries for TODO and the records of the ledger are left to the landing.

What moved, and why. glue_toggle of the BCF baseline, which holds the treatment forest's
coefficients: its forest's sd is now exact where the map divided it by 0.674. Every other recorded
scenario is bitwise: 55 of 55 against equivalence-32418b6e under `--strict-coverage`, 14 of 15
against bcf-equivalence-1b7d730c, 11 of 11 against multinomial-equivalence-80b1c8d4, and the four
snapshot files unchanged, on a `--preclean` reference install of the final code. bartCause's `bcf()`
keeps its draws at its default and with both sds stated, and moves with `update.a = FALSE` or
`update.b = FALSE`.

The pair script, 42 rows of the 29 above, old side on 32418b6e's code and new side on push 1's: 0
unexpected. Rows 23 and 24 are not only equal to 1e-10 but identical; row 30 is a constant response,
whose copy the base build refuses.

Gates on the tip: tests/cpp 401 lines ok, and under ASAN and UBSAN 401 with 0 diagnostics; the full
suite through `test_package(at_home = TRUE)` 27454 results, 0 failed; the 43 exact-gate commands of
exact-gates.yaml in quick, all passing; the composition matrix 269 checks, 0 disagreements; lintr,
air, rc-codoc, win-drift, doc-freshness and verify-anchors clean; dbarts.h no diff; `R CMD check
--as-cran --no-manual` on a clean tarball one NOTE (the Date field).

Size: 2016 lines added and 871 removed over 43 files after the recheck's commit, the baseline file
aside: src 187 and 90, R 80 and 102, tests/cpp 380 and 69, tinytest 660 and 323, man 75 and 87,
benchmarks 406 and 33, docs and TODO 223 and 163. Package code 459 lines against the stop of about
900, the whole 2887 against about 5000. bartCause 41 and 15.

Mutants, each applied alone to a clean export of the tip, installed with `--preclean`, and run against
tests/cpp where it is the engine's, the named tinytest files and, where marked, forest-prior-exact.R in
quick. Every one is killed:

| mutant | killed by |
|---|---|
| the kind not carried from the data door (every basis numbers) | test-forest-kind.R and test-forest-held.R stop: a held factor is refused as a held number |
| the kind lost on a copy (the copy's levels dropped) | test-forest-kind.R, 3 |
| a held value left by position | tests/cpp 6; test-forest-held.R 46 |
| a held forest with no basis second held at 0 | tests/cpp 3; test-forest-held.R 8 |
| the 0.674 left on a held factor | tests/cpp 2; test-forest-held.R 21; the gate, 20 functionals |
| a held forest with no basis left at the anchor whatever its sd | tests/cpp 3; test-forest-held.R 16; the gate, 20 |
| a held forest with no basis given one unit by default under gaussian | tests/cpp 2; test-forest-held.R 2; the gate, 12 |
| a leaf's sd not divided by the root of the tree count | tests/cpp 91; test-forest-held.R 39; the gate, every arm |
| sqrt(2 / K) left in the default of a forest with a basis | tests/cpp 56; sd unit 35; bcf family 16; the gate, 5 in the arm of three forests |
| a held forest with no basis at 2 sd(y) by default | tests/cpp 2; test-forest-held.R 2; the gate, 18 |
| `amplitude.prior.variance` taken and ignored | test-forest-held.R 3 |
| a write to a held forest with no basis that leaves its leaf scale | tests/cpp 1 |
| a held block installed from a state | tests/cpp 1; test-state-not-model.R 1 |
| an interim refusal of a held forest with no basis second left standing | test-forest-held.R and test-forest-defaults.R stop |
| the anchor of a constant response left at 0 | tests/cpp 3; the sd unit file stops |
| a held factor at (1, 0) | tests/cpp 6; test-forest-held.R 66 |
| the bridge reporting every forest drawn | test-forest-held.R 59 |
| the reader's " (prior median)" dropped | test-calibration-midchain.R 2 |
| no throw for a forest of levels or numbers handed no basis | tests/cpp 1 |
| no throw for a forest with no basis handed one | tests/cpp 1 |
| the held default at the node scale, not over k (half the range) | tests/cpp 2; test-forest-held.R 2; the gate, 18 |
| a held forest with no basis given a half-Cauchy median | tests/cpp 1 |

Not run: the variant that takes a quarter of every row's range where rows are out (test-forest-held.R's
equality with `bart()` at ten rows of weight zero far out stands for it, the gate having no weights).

Fix round, from the review (one must-fix, four should-fixes). The held rule's exact gate is renamed
forest-prior-exact.R, and its usage line, its printf, a test comment and forest's help say prior. Three mutants the review found surviving
are killed by new tests, each run on a clean export with `--preclean`:

| mutant | killed by |
|---|---|
| the held default taken from a fresh anchor, not the retained one | test-forest-held.R 2: a kept-scale `setResponse` then `copy()` |
| a hand-written bcf-shaped state installing a held b | tests/cpp 1: `restoreGlue` of a state with no widths, b held |
| a write to a held factor resolved as drawn | test-forest-held.R 2: `setLeafPrior` of an sd on a held factor |

The sd unit note says push 1 made a forest with a basis one unit whatever the number of forests; the
BCF baseline's MANIFEST row names its oracle's metrics, a largest elementwise relative difference of
9e-14 and a mean of at most 2.2e-15; bartCause's held default reads the family (0057d2a, with a test), and a test pins its weighted rows (4a70e7d).
The branch was rebased onto 33f9c299, which moved no draw of the three baselines, and the BCF baseline
renamed to the rebased commit.

## Landing note

Push 1 landed on 2026-10-17 as f9bc3636..5a62933e, fast-forwarded onto 33f9c299 after a second rebase
that changed only cites, with bartCause's 30cd6fc, 0057d2a and 4a70e7d fast-forwarded onto dbarts-1.0.
The calls made planning and building it are dec-A216. The x86-64 leg ran on the code landed: tests/cpp
402 ok plain and under ASAN, the suite 27617 results with 0 failed, the equivalence compare 55 of 55
statistically indistinguishable, and the bcf and multinomial compares differing bitwise only where the
build before does on x86. Push 2 must land before 1.0-0, which removes `amplitude.prior.scale` and with
it the prefix match of `$amplitude`.

## Push 2 as built

Built 2026-10-17 on bartcore at 738b1015 (push 1 landed), not landed. Commits, in order: the relative
bound, step 2.1, the first commit (a51ac7fe); the engine and the bridge, steps 2.2 to 2.4 with their
tests/cpp (f2ccd656); R, steps 2.5 to 2.7 and 2.10, with test-forest-scale-record.R,
test-forest-leaf-prior-reader.R and the units twin (c0b58aaa); forest-prior-exact.R's numeric arms, 2.8
(fcbba19b); the help, 2.11 (1dfa5989, the last commit that changes package code); the four BCF scenarios,
2.9 (8b48cccf); the design note and records, 2.11 (f1fd0de7); a pin of a numeric coefficient's start
standing second (3e241cfc); the BCF baseline (0f942fad). No consumer is edited: bartCause reads no entry
this push removes, and stan4bart, treatSens and bairrtt declare one forest.

Against the steps:

- 2.1 as planned, the comparison `<=` so that a forest whose multipliers are all zero loses every row.
  tests/cpp `testRelativeRowBound`; the tinytest twin is test-forest-units-twin.R's first arm. The
  seam fixture's band edge moved from 2^-27 to 2^-29, its b1 being 0.25. Gated alone: tests/cpp 403 ok;
  the suite 27538 results, 0 failed; on a reference build the four snapshot files passing and 55, 15 and
  11 scenarios identical. Against its parent,
  `bench-sampler.R record` on 738b1015 then `compare` on a51ac7fe, run by the orchestrator on the idle x86
  bench box (4 cores, load about 1), one machine: no metric regressed more than 5%. The ratios: embedded
  offset gibbs step 1.000, binary run 0.988, run t200 1.000, run t75 1.000, run n10000 0.994,
  setPredictor accept 1.035, setPredictor reject 1.000.
- 2.2 as planned: [`basisColumnScale`](../../src/bartcore/combiner.hpp) and
  [`BasisScaleError`](../../src/bartcore/combiner.hpp); the chain's rows in `rowsForBasisScale`, the
  status in `ForestCalibration::basisScaleFromAllRows`. tests/cpp `testBasisColumnScale`.
- 2.3 as planned: `ForestSpec::basisScale`, the numeric row of `forestMapEntry` (variance 1, start 1,
  held at 1), `basisRowNorm` gone, the constructor calling `mapLeafScale`, a swap deriving nothing.
  `ForestCalibration` no longer carries a factor, a divisor or a row norm; tests/cpp reads the entry
  through `TestPeer::mapEntry`. tests/cpp `testNumericMultiplierLaw`.
- 2.4 as planned, and `bartcore_create` takes a fifth argument, the mask (below).
- 2.5 as planned: the slot, `dataBasisScale`, `recordBasisScales` in the sampler's `initialize`, the
  entry dropped in `resolveSamplerSpec` where a declaration replaces a basis.
- 2.6 as planned; the refusal is worded by the bridge, which reads the forest's label from the forests'
  record.
- 2.7 as planned; `prior.sd.of` for a numeric forest is "forest, per unit of basis", with " (prior
  median)" where drawn.
- 2.8: nine arms of one column (a 0/1 column drawn and held at its default and drawn at 0.6 of the unit;
  a column of four values, sd about 0.11, drawn and held at the default and at 3 of the unit; a fifth of
  the rows at weight zero, their column at 2, drawn and held) and two arms of three trees, the number
  drawn and held. The coefficient's quadrature is a grid on [-8, 8] at 0.005 under its prior.
- 2.9 as planned: numeric_default, numeric_swap_write_copy, numeric_held and held_plain, with literal
  seeds.
- 2.10: the respell touched 15 test files. A creation on a column constant over the rows in states an
  sd (test-bcf-family.R's first forest of ones, through its `statedOnFirst`, and
  test-response-scale-rows.R's z, which is 1 on every row in), or is turned into a refusal
  (test-forest-kind.R's constant column, at the seven doors).
- 2.11 as planned.

Calls made building, where the plan was silent or wrong:

- Under probit and logistic the 0/1 weights' mask was installed after the sampler was created, where
  the engine could not read it. The mask is handed to creation, `new("dbartsSampler", control, model,
  data, active = )`, which installs it as `dbarts()` did and passes it to the bridge, where
  `AmplitudeSpec::activeRows` gives the rows of a column's scale and nothing else reads it. A mask
  installed after creation moves no prior, as rule 3 has it; dbartsSpec's help says how to build the
  sampler `dbarts()` builds. Rejected: taking the scale again at the first mask install, a scale taken
  twice; and every row under the latent families, against settled call 1.
- One warning names the response and the column: `dbarts()` holds the spec's response warning and
  raises it only where the engine reports no column read from every row, in which case the sampler's
  `initialize` raises the joint one. Rejected: two warnings; R deciding the column's fallback.
- The bridge keeps its raw `amplitude.prior.variance` and `amplitude.prior.scale` columns, which the
  reader does not report, so a test can read a drawn forest with no basis's half-Cauchy median.
- The engine refuses a swap handing a numeric forest other than one column, and the bridge raises it.
- `summary()`'s line takes each forest's number unnamed, so a numeric forest's column name does not
  make its label `dose.dose`.
- The BCF baseline is named after the last commit that changes package code, and its meta.rev set to
  it: 1dfa5989 as built, 84c08a20 after the fix round, 20eed7ed after the second.
- A by-hand probit sampler of the plan's 2.5 "which rows" is built with the mask at creation.

What moved, and why. Every model with a numeric multiplier, as the `rng:` line says; no recorded
scenario has one. On a `--preclean` reference install of the final package code: 55 of 55 against
equivalence-32418b6e under `--strict-coverage`, 15 of 15 against bcf-equivalence-47838362, 11 of 11
against multinomial-equivalence-80b1c8d4, the four snapshot files unchanged; then 19 of 19 against
bcf-equivalence-1dfa5989 as recorded.

The pair script, 45 rows, old side on 738b1015 and new side on push 2: 0 unexpected. Rows with a
numeric column differ (03, 15, 17, 19, 26a, 29 with a number, 32); rows of factors and forests with no
basis, held or drawn, are identical (29 with factors only included); a held numeric column is created
(27); a constant column at its default is refused (31); dose and 2^-27 dose draw the same on the new
side and not on the old (28).

Gates on the tip: tests/cpp 405 lines ok, and under ASAN and UBSAN 405 with 0 diagnostics; the full
suite through `test_package(at_home = TRUE)` 28001 results in 261 files, 0 failed; the 43 exact-gate
commands of exact-gates.yaml in quick, all passing, forest-prior-exact.R's largest |z| 3.17; lintr,
air, rc-codoc, win-drift, doc-freshness and verify-anchors clean; dbarts.h no diff; `R CMD check
--as-cran --no-manual` on a clean tarball one NOTE (the Date field); the R-loaded address sanitizer
over the eight test files the change reaches, 0 diagnostics. bartCause's suite on dbarts-1.0 at
4a70e7d, unedited, against private installs of both packages: 1367 expectations, 7 failing in the 5
tests that fail the same way against the base build (three regression comparisons with old results, a
grouped predict and a stan4bart predict). stan4bart, treatSens and bairrtt read nothing this push changes
(searched, not run).

Size: 3526 lines added and 1110 removed over 57 files, the baseline file aside: src 377 and 157, R 232
and 99, tests/cpp 752 and 174, tinytest 1353 and 432, man 123 and 97, benchmarks 430 and 8, docs, TODO
and workflows 259 and 143. Package code 865 lines against the stop of about 1100, the whole 4636 against
about 6400.

Mutants, each applied alone to a clean export of the tip, installed with `--preclean`, and run against
tests/cpp where it is the engine's, the named tinytest files and, where marked, forest-prior-exact.R in
quick. Every one is killed:

| mutant | killed by |
|---|---|
| the bound left absolute | tests/cpp 6; the units twin 4 |
| the bound relative to the largest multiplier of every forest | tests/cpp 5 |
| a zero multiplier kept in the update | tests/cpp aborts on the cache check |
| the default per unit (the scale read as 1) | tests/cpp 37; the record 3; the reader 2; the gate 80 |
| a scale 3 percent off | tests/cpp 23; the record 17; the gate 5 |
| n for n - 1 | tests/cpp 22; the record 17 (the gate does not see it) |
| the scale over the nonzero rows | tests/cpp aborts on a 0/1 column |
| the scale over every row where rows are out | tests/cpp 1; the record 12; the gate 16 |
| one pass for the mean | tests/cpp 8; the record 1 |
| constant decided by the computed sd | tests/cpp 8; the record 20 |
| the near-constant bound dropped | tests/cpp 1; the record 1 |
| a numeric coefficient's variance left at 1/2 | tests/cpp 15; the gate 45 |
| a stated sd divided by the scale | tests/cpp 11; the reader 4 |
| the scale taken again at a re-creation | tests/cpp 1; the record 3 |
| the record kept where a declaration replaces the basis | the record 1 |
| an unread record checked | the record stops |
| a read record not checked | the record 2 |
| a numeric coefficient started at 0 as the second forest | tests/cpp 1 |
| a held numeric column at 0 | tests/cpp 3; test-forest-held.R 63; test-forest-arguments.R 4; the gate 47 |
| a swap that takes the scale again | tests/cpp 2; the record 1; test-forest-basis-r5.R 1 |
| a numeric forest's entry left unnamed | the reader 9; the record 2 |
| the prior-predictive refusal dropped | the record 1 |
| the mask not handed to creation under probit | the record 6 |
| the fallback to every row dropped under probit | tests/cpp 1; the record stops |

Fix round, from the review (no must-fix, three should-fixes, seven nits taken). Commits: the tests
(536b9d60), the help (84c08a20, now the last commit that changes package code), the design note and
TODO (9e1f9bb3), the BCF scenario (7df99a10) and its baseline (5e90c526).

- numeric_swap_write_copy wrote the sd before it copied, so the copy was a stated forest and never read
  the record it was handed: clearing `data@basis.scale` before `copy()` left its draws identical. The
  copy is now made after the swap and run at its default before the write. The same probe now moves
  the draws (the doubled sd 90.1 against 9.0), and the harness, with the record cleared in that
  scenario, fails it bitwise and at a largest |z| of 11.9. Re-recorded as bcf-equivalence-84c08a20:
  against 1dfa5989's file only numeric_swap_write_copy moves; the 15 earlier scenarios are 15 of 15
  identical against 47838362, restricted's varcount.tau differing in its dimnames attribute alone.
  1dfa5989's file is dropped from the branch.
- The units twin gains a column of size about 3 with one row at 0.9 times 2^-26 of its largest, the
  column against 2^-20 times it, so the bound is tested above 1 as well as below.
- test-forest-scale-record.R gains a column holding 1e9 on its rows of weight zero, whose scale is
  that over the rows in, and the gaussian joint warning's text in full.
- tests/cpp: a numeric forest refuses a swap of two columns in the engine.
- The help: `new("dbartsSampler", ..., active = )` on the class page, read at creation under probit
  and logistic and installed afterwards under every family; the two warnings of the by-hand route under
  gaussian, in dbartsSpec's `active`; and in forest's `sd`, probit and logistic reading no sd(y).
- The design note's consumer case is gone, along with its TODO entry: the shipped header has no create
  entry and cannot describe a sampler of several forests, so the case cannot arise. The TODO entry the
  plan lists for it is not added.

The review's four surviving mutants, each run on a clean export with `--preclean`:

| mutant | killed by |
|---|---|
| the bound `tol * min(largest, 1)` | the units twin 1 |
| the near-constant check's largest read over every row | the record stops (refused as constant to rounding) |
| the engine's numeric swap width check dropped | tests/cpp 1 |
| the gaussian joint warning dropping the column | the record 1 |

Gates on the tip: tests/cpp 405 lines ok, and under ASAN and UBSAN 405 with 0 diagnostics; the full
suite through `test_package(at_home = TRUE)` 28004 results in 261 files, 0 failed; the 43 exact-gate
commands in quick, all passing; lintr (0 lints), air, rc-codoc, win-drift, doc-freshness and
verify-anchors clean; on a `--preclean` reference install 55 of 55 against equivalence-32418b6e under
`--strict-coverage`, 11 of 11 against multinomial-equivalence-80b1c8d4, 15 of 15 against
bcf-equivalence-47838362 and 19 of 19 against bcf-equivalence-84c08a20, and the four snapshot files
passing (14, 3, 7 and 3); `R CMD check --as-cran --no-manual` on a clean tarball one NOTE (the Date
field). Size of the round: 65 lines added and 27 removed over 16 files, the baseline files aside.

Second fix round, on the reviewer's E7, ruled by the orchestrator: a row out of the likelihood must not
move the model. `formForestResponse` took the forest's largest multiplier over every row, so a value
of 1e9 on rows of weight zero put every row of positive weight under 2^-26 of it, and the forest fit
from its prior alone. The largest is now taken over the rows of positive weight, every row where none
has one, and a forest whose largest is zero loses every row (20eed7ed). A row of weight zero is
treated as before in the second loop.

- tests/cpp `testRelativeRowBound`: a row of weight zero at 1e9 moves no row of positive weight out; a
  column of zeros on every row of positive weight loses every row; with every weight zero the rows
  leave as with no weights.
- test-forest-scale-record.R: the column with 1e9 on its ten rows of weight zero draws what the same
  column with 5 there draws, in sigma, the variable counts, the coefficients and the fits at the rows
  in. A fit at a row out is read back through a residual its multiplier's size rounds, so it differs
  there in the seventh digit.
- The mutant, the largest over every row: tests/cpp 2; the record 1, its sigma, counts, fits and
  coefficients all moved.
- No recorded draw moved: 55 of 55, 11 of 11, and 15 of 15 against 47838362 and 19 of 19 against the
  84c08a20 file, so the BCF baseline is renamed bcf-equivalence-20eed7ed with its meta.rev set, not
  re-recorded (8e7a33fd).
- The step 2.1 compare covered the bound before this change; the orchestrator reruns `bench-sampler.R`
  for the touched hot path on the bench box.

Gates on the tip: tests/cpp 405 lines ok, and under ASAN and UBSAN 405 with 0 diagnostics; the full
suite through `test_package(at_home = TRUE)` 28005 results in 261 files, 0 failed; the 43 exact-gate
commands in quick, all passing; lintr (0 lints), air, rc-codoc, win-drift, doc-freshness and
verify-anchors clean; on a `--preclean` reference install the four snapshot files passing (14, 3, 7
and 3) and 55, 11, 15 and 19 scenarios identical, the last against bcf-equivalence-20eed7ed as
renamed. Size of the round: 84 lines added and 8 removed over 8 files, this record aside.

## Landing note, push 2

Push 2 landed 2026-10-17 as a51ac7fe..5a1e9dca on 46e1a644, 20 commits: the build (a51ac7fe to 0f942fad), its
record (96fcde77), the fix round (536b9d60 to 486a7893), the second fix round on E7 (20eed7ed to 0763eed8),
and at the landing the BCF baseline renamed to bcf-equivalence-20eed7ed, the last commit that changes package
code once rebased, with the hashes this plan cites (5a1e9dca). It was built on 738b1015 and rebased onto
46e1a644, the records of warm-start slice C, without conflict. The calls are dec-A219.

Gates on the rebased tip: the suite through `test_package(at_home = TRUE)` 28159 results, 0 failed; lintr (0
lints), air, rc-codoc, win-drift, doc-freshness and verify-anchors clean; the BCF compare 19 of 19 bit for bit
against bcf-equivalence-20eed7ed. The bench compare of the second fix round on the x86 bench box, 1ef44d1e's
parent against it before the rebase: no metric more than 5% slower (run t75 1.016, setPredictor reject
1.012, the rest within 1%). The x86-64 leg on the second fix round's tip: tests/cpp 405 ok, plain and under
ASAN with the leak sanitizer; the suite 28038 results, 0 failed; the equivalence compare 55 of 55
statistically indistinguishable; the BCF and multinomial compares within the cross-host tier 1 bound, 19
of 19 and 11 of 11; and against a BCF baseline recorded on the box from 46e1a644, 15 of 15 bit for bit.

Not run: Windows; bartCause's suite after the rebase (run on the build, 7 failing as on the base build).

Left: push 3, the printed block.

## Push 3 as built

Built 2026-10-17 on bartcore at 5a1e9dca (push 2 landed), rebased onto a000d239 (its records), R only, not landed. Commits: the block, step
3.1, with test-forest-print.R, step 3.2 (5625cd90); the help and the design note, step 3.3 (2510b0fa).
Nothing under src/ and no draw moves; `dbarts.h` has no diff. The block described below is the full one;
[Push 3 fix round](#push-3-fix-round) trims what `print` and `show` give to a short block.

Against the steps:

- 3.1 as planned: [`forestPriorLines`](../../R/forestPrint.R) builds the block from the reader's entries
  and, for a numeric forest, the column as the object holds it. [`fitSynopsis`](../../R/generics.R)
  prints it last, after the monotone prior where there is one, so the tree counts, which a fit with no
  kept sampler lacks, are where they were; `show` prints it after the call; a sampler created with
  `verbose = TRUE` prints it once, after the engine's summary.
- The plan said nothing more is stored on a fit. The header's number is not among the reader's entries
  and is not the response's sd where the response is constant or rows are out, and R recomputing the unit
  would be the copy the sd unit slice forbids. The fit therefore keeps one number: the unit the engine
  recorded on the sampler's control, `attr(fit, "forest.unit")`, beside `forest.labels` and, like it,
  not an element of the list. Rejected: reading the unit off a kept sampler only (a default fit keeps
  none); sd(y) of the stored response (wrong on a constant response and on rows out).
- The rows in the likelihood of a fit are its stored `weights` and `active`; of a sampler, its weights
  and `activeRows`.
- 3.2: test-forest-print.R, 36 results. The block by exact string for the model of four forests at
  `show`, creation and `print` of a fit with and without a kept sampler; a stated forest, a held one, a
  held forest with no basis; probit and logistic headers; a swap to dose / 10; rows of weight zero
  holding 1000, at creation, on a fit and after `setWeights`; `scale()` printing mean 0; one value a
  thousand times the rest in the range; a constant column with a stated sd printing no "per sd"; every
  number of the numeric lines read back from the text against the reader and R's `mean`, `sd` and
  `range` to four digits; a copy printing nothing; one forest printing nothing; no warning in the file.
- 3.3 is in man/forest.Rd's `sd` item, with the two things to read, and man/bartBT.Rd, where `print.bart`
  is described (the plan named man/bart.Rd, which does not describe it).

Calls made building:

- Numbers print as `format(x, digits = 4)` gives them: at least four significant digits and every digit
  of an integer part, so 58203 and not 5.82e+04. The plan's figures (`0.0827`) are hand-rounded and not
  a rule. Rejected: three digits; `signif` then `as.character`.
- A column is constant in the block where its sd is under 1e-11 of its largest value, the engine's
  bound, not where it is exactly 0.
- A "default" line ends with the recorded sd where the column's sd in force differs from it by more than
  1e-8 relative; a stated line never does.
- A copy does not print: `copy()` holds the block off, so a verbose sampler prints once per creation.
- The verbose creation prints the block from R, after the engine's own summary, which is not edited:
  the plan puts nothing under src/.

Mutants, each applied alone to a copy of the tree, installed, and run against test-forest-print.R. Every
one is killed:

| mutant | killed by |
|---|---|
| a column's printed sd read from the record | 2 |
| the per-sd figure from the recorded sd | 2 |
| the range left out | 18 |
| a label by position | 17 |
| a mean not passed through `zapsmall` | 2 |
| every row read, not the rows in the likelihood | 1 |
| "per sd" printed for a constant column | 1 |
| a copy printing the block | 1 |

Gates on the tip, which is on a000d239, after a `--preclean` install: the full suite through
`test_package("dbarts")` 28100 results in 262 files, 0 failed; lintr (0 lints), air, rc-codoc, win-drift
and doc-freshness clean; the 43 exact-gate commands of exact-gates.yaml in quick, all passing, run
because a fit now carries the unit. The bitwise compares, the reference build and the sanitizers are not
run: no file under src/ changes and no draw moves, which the R-surface tier does not ask them for.
Size: 635 lines added and 5 removed over 8 files before this record: R 214 (of which the block 201),
tinytest 391, help 15, design note 15. Package code 214 lines against the plan's 520 in all.

### Push 3 fix round

Built 2026-10-22 on 3a9559db, R only, from one opus review (ae70169a) and the package owner's ruling of the
same day on the printed block. Nothing under src/ and no draw moves; `dbarts.h` has no diff. Commits: the
code and tests (2ef87643); the help (86f7a908); this record.

What changed, by the review's findings:

- The gaussian header reads `forests, sd in the response's units (sd(y) 1.452):`, the help's own symbol,
  whose `sd` item covers rows out of the likelihood, the offset, the fallback to every row and a constant
  response. The label built, "its sd", was false on an offset fit (the number is sd(y - offset)), on a
  constant response (a quarter of the window of width one) and after `setResponse` (the creation's number,
  against the response in force). Latent families keep `(error sd 1)` and `(error sd 1.814)`. Pinned on an
  offset sampler and fit, a constant response (the creation's one warning, none at show, print or summary)
  and a response replaced by `setResponse(10 * y)`, which leaves the number and the defaults the creation's.
- A fit made with `keepFits = FALSE` keeps no bases, labels or unit and gives no block; the help says so in
  `forest` and `bartBT`, and a test pins the absent elements, the absent block and, in its summary, the sds
  kept on the line under the table. Rejected: a block without the column figures, which is a second format
  for a fit holding no column.
- The per-sd figure is the size of the forest's term over one sd of the column (`sd` times the column's sd),
  not the size of one sd of the column; the help said the second.
- The reading of a fit with no unit is gone: `forestPriorLines` takes `unit` with no default, the `bare`
  test and the record's "a fit without it (made before)" go. Only shipped versions are owed compatibility.
- `attr(fit, "forest.unit")` is documented beside `forest.labels` in the Value sections of bart.Rd and
  bartBT.Rd. Whether a fit carries it at all went to the package owner, who ruled it stays.
- Help wording: "a mean far from 0 is a column that should be centered" now holds for a dose too (a column
  that never comes near 0 should be; a dose, where 0 means none, is left); "a default taken from a column
  that has since been replaced says so" names the weights, the active rows and the fallback to every row.
- Numbers print as `format(x, digits = 4)` gives them, which the help and the record now say; the test
  helper `g()` is that call, where `formatC(.., "g")` differed from 1e4 up (5.82e+04 against 58203).
- `show`'s call line printed `dbarts(formula = y ~ x1 + x2, ...)` for a model written with `forest()`
  terms, since `dbarts()` replaced the formula by the one its terms were taken out of before it recorded the
  call. The sampler now records the formula as written (8 lines in [`dbarts`](../../R/dbarts.R)), pinned for
  a formula written out and one passed by name.
- Tests: test-forest-print.R 36 results to 132, and test-fit-stores-k.R's summary line, which named each
  forest's sd, now reads the block. The `zero` and `single` verbose creations are captured, the output and
  the messages, so the run prints nothing. Counted warnings wrap every show, print and summary of a fit as
  well as every creation. New pins: probit and logistic on a fit; under probit 0/1 weights are the mask, a
  value of 1000 on a row of weight zero stays out, on a fit and on a sampler; a stated line never takes the
  "taken from" mark after a swap then `setLeafPrior(forest(sd = 30))`; a column of 61.7 + 1e-13 sin(i)
  reads as constant; the hold released after a copy; the index padded at ten forests; the fallback to every
  row under probit with every weight 0, with its warning; numbers from 1e4 up.

The owner's ruling of 2026-10-22, the trim. `print` of a fit and `show` of a sampler give a short block: one
line for each forest, its label, its basis, for a numeric basis the column's mean and sd over the rows in
the likelihood, and the sd it resolved to (per unit, or between the two levels of a factor), with no range,
no per-sd figure, no source or "taken from" mark and no mark of a held coefficient. Its header is plain,
`forests:`, and its last line points to the full block. `summary()` of a fit gives the full block as built,
header and all, after its table; a verbose creation gives it too; a copy prints nothing; a fit of one forest
prints and summarizes as before. The line under summary's table leaves out the `leaf.prior.sd` entries where
the block gives them and keeps any others.

Calls made building it, with the alternatives:

- Header `forests:` with the last line `(each forest's prior in full: summary())`, the owner's example.
  Rejected: a header that says "short"; "see summary() for more".
- A sampler has no `summary`, so its last line reads `(each forest's prior in full: $getLeafPrior())`, which
  returns each forest's entry. Rejected: a `summary` method on the sampler, which is a new public surface;
  one footer naming `summary()` for both, which would point a sampler's reader at nothing.
- In summary the block follows the table and its notes (the line under it, the starting-sigma line, the
  R-hat note) after a blank line. Rejected: between the line under the table and the R-hat note, which
  separates the note from its table; before the table.
- The summary object carries the block as `forests`, absent where there is none; its `fixed` still holds
  `leaf.prior.sd` and extract is untouched, and only the printed line drops it. Rejected: dropping it from
  the object, which would move data to move text. A fit with no block (`keepFits = FALSE`) keeps the entries
  on the line, as there is nothing else that gives them.
- Not done, and the owner's to rule: after the trim nothing prints the full block of a live sampler after
  its creation, so the "taken from a column of sd" mark can only show at a creation from a data object whose
  column was replaced. The tests reach the full block of a changed sampler through
  [`samplerForestLines`](../../R/forestPrint.R). Alternatives: show the mark in the short line where it
  applies; drop the mark and its code.

Mutants, each applied alone to a copy of the tree, installed, and run against test-forest-print.R, the
tests that fail named by result number (the file now has 132). Every one is killed. The review's:

| mutant | killed by |
|---|---|
| A: a fit's family passed as "gaussian" | 49, 53 (the latent header at summary of a probit and a logistic fit), 114, 120 |
| B: a fit's mask left out | 113, 114 (the 1000 on rows of weight zero, at print and summary of a probit fit) |
| C: a stated line takes the mark | 61 (a swap, then `setLeafPrior(forest(sd = 30))`) |
| D: constant only where the sd is exactly 0 | 80 (the column 61.7 + 1e-13 sin(i)) |
| E: no `forest.unit` on a fit | 13, 14, 18 (the attribute, the full block at summary) and 49, 53, 66, 99, 103, 114, 120 |
| F: the header's unit from sd of a fit's response | 49, 53, 66, 99, 103, 114, 120 |
| G: a sampler's mask left out | 109, 110 (a probit sampler, at creation and show) |
| J: the hold not released | 35, 36 (the creation after a copy) and 26 others |
| K: the index not padded | 89 to 91 and 93 to 95 (ten forests, full and short) |
| L: no fallback to every row | 117 to 120 |
| M: a warning in a fit's block, or a sampler's, or a summary's | 132 (the count at the end), all three |

Added for this round's changes:

| mutant | killed by |
|---|---|
| the call recorded with the stripped formula | 4, 5, 7 |
| the header read "its sd" | 13 tests, among them 8, 14, 63, 98, 102 and 105 |
| the header's unit from sd of a sampler's response now | 9 tests: 46, 50, 63, 98 (offset), 102 (constant), 106 (the replaced response) and 66, 109, 117 |
| `formatC(.., "g")` trimmed, as the old helper was | 82 to 86 (the column of 1000 ages) |
| a fit's weights left out | 65, 66 |
| the short block with the range, the per-sd figure, the source or the held mark | 29, 24, 25 and 6 tests |
| the short header or footer wrong, the fit's pointer or the sampler's changed | 36, 18, 8 and 9 tests |
| the summary without its block, or with the sds still on the line | 10 tests; 16, 17, 22, 23 |
| a verbose creation short, a show or a print full | 33, 28 and 10 tests |
| the summary's block short | 12 tests |

Gates on the code tip (86f7a908, on a000d239) after a `--preclean` install, mode "shipped", the install newer
than the sources: test-forest-print.R 132 results and test-fit-stores-k.R 500, 0 failed; the full suite
through `test_package("dbarts")` in one process 28199 results in 262 files, 0 failed (28100 before: 96 more
in test-forest-print.R, 3 in test-fit-stores-k.R); the 43 exact-gate commands of exact-gates.yaml in quick,
39 of the main loop and 4 of the monotone job, all passing, run because a fit now carries the unit and a
summary a new element (the two cross-host compares not run: no draw moves); lintr (0 lints), air, rc-codoc
(48 methods), win-drift and doc-freshness clean, each on its own exit status; `R CMD check --as-cran` on
a tarball built from a clean `git archive` of 86f7a908 staged outside the tree, 1 NOTE ("Date field over a
month old"), its tests OK; on a reference install (`mode` "reference") the four snapshot files, binary 7,
multithreaded 3, single-threaded 14 and xbart 3 assertions, 0 failed. This record is the only change after
that tip. Size over 3a9559db: R 157 lines changed over three files, tinytest 688, help 56.
