# hurdle-completion: a test set, weights, subset, offsets and the formula interface on hurdle.lognormal fits

Status: PLANNED 2026-11-03 on main (170635476). Not built; a blind critique follows. Two open calls.
Rulings, quoted where they settle a call: dec-A15 (fitted takes ci.level after sample; of the hurdle's missing test channel,
"make sure that its missing test channel is treated as part of feature-completion work and we don't forget it"), dec-A114
(hurdle predict's offset is a reserved slot; "Hurdle offsets are coming."), dec-B07 (one predict argument order), dec-B236
(fitted takes what extract takes), dec-A19 and dec-B82 (one token, composed in R from two independent samplers), dec-B14
(the hurdle's sampler-only mode stays refused), dec-A105 (a hurdle refuses a named sd, "whose two halves would read one
number in two units"), dec-A112 and dec-A113 (missing values at a fit-time test set and on zero rows only), dec-B154,
dec-B184, dec-B340 and dec-B341 (offsets: a term and the argument alike, re-evaluated on newdata and added to there,
replacing the fit's at the training rows), dec-A111 (offset.test defaults to the training offset where the rows match),
dec-A137 and dec-B349 (an own-class predict takes an offset as every other family's does), dec-B351 (each argument that holds
rows is evaluated once per fit), dec-B97 (the survival interface finished in R, "A user fits a survival model the way they fit
any other"), dec-B563 (frequency weights; probit, ordinal and aft out of scope; "hurdle's weights stay with its feature
completion (dec-A15)"), dec-B13, dec-B564, dec-B565, dec-B569, dec-B570, dec-B580 (a 0 row's predictive draws are NA on
every family; "0 weight means out of sample"). The runs behind it are the probes under scratch/hrd-plan/probes, on a shipped
build of main 8fce28568 (scratch/libs/hrd-plan), 0.9-34 (scratch/libs/cran) and pscl 1.5.9 and countreg 0.3.0
(scratch/hrd-plan/landscape-lib); brms 2.23.0 was read from its installed source.

agent: S1 to S4 one sonnet implementer each (effort medium); S5 sonnet. One sonnet review for S1 and S2 together and one for
S3 and S4, each reviewer checking its fix round; S5 takes the R-surface tier's review as a gate-only slice
([Process by risk](README.md#process-by-risk)). Tier "R surface, no draw change": no engine, bridge, header or state-format
file moves, and no call taken today changes a draw.
rng: neutral. Every hurdle call taken today (the matrix interface with none of test, weights, subset, offset) keeps its draws
bit for bit: the equivalence harness's hurdle scenario, hurdle-reduction.R, hurdle-exact.R and test-hurdle.R's pins do not
move. New calls take new draws, each part held draw for draw to the standalone probit or gaussian fit it is.
window: before 1.0-0 (dec-A15). One branch, hurdle-completion, merged once (dec-B494). R only, so it need not wait on
coefficient-prior-live or the engine queue; order against frequency-weights below.
budget: about 2,200 lines: R ~690, help ~140, NEWS ~5, tinytest ~1,150, gate scripts ~200, records ~40. Agent time,
estimated and not measured: S1 ~6, S2 ~6, S3 ~5, S4 ~8, S5 ~5, S6 ~1, three reviews with fix rounds ~11, gates and check ~3:
~45. Plans here run 1.5 to 2 times low, so 70 to 90.

## Summary for the maintainer

A hurdle fit today takes only a matrix of predictors and a response. It refuses a formula, a test set, weights, a subset
and an offset, and its fitted method has no way to ask for test rows. It is the only family that refuses a formula. This
plan finishes it so that a hurdle fit is made the way any other fit is made.

What a user gets. bart(y ~ ., data, family = "hurdle.lognormal") works, with subset, weights, an offset and a test set
written as on any other family. The test set gives stored draws at those rows, so extract and fitted take sample = "test"
and a fit kept without its trees still answers for them. Both parts always see the same rows: a subset or a missing value
removes a row from both, as one model frame does in lm and pscl's hurdle. predict takes an offset in the place held for it
and, as a gaussian or probit fit's predict does, a weight of 1 to draw a held-out row as a new observation.

Two choices are yours.

First, which part an offset enters. The zero part is a probit and the positive part models log y, so one number cannot
mean the same thing in both, the reason a single leaf-prior spread is refused on this family. pscl and countreg put their
offset argument in the count part, and brms puts an offset in the part whose formula carries it, the size part for the
main formula. The package already reports the positive part as a hurdle fit's link scale. I recommend that a plain offset
enters the positive part, and that a list naming zero and positive gives each part its own, using the part names you
chose for this family's k and variable counts.

Second, what weights it takes. pscl and brms read a hurdle's weights as counts of identical rows in both parts. The zero
part is a probit fit, and a probit takes only 0 and 1, which say whether a row is in the data; counts on probit were
scheduled after 1.0-0 in August and left unscheduled when you ruled frequency weights. I recommend the hurdle takes what
its probit part takes: 0 and 1 now, so a row can be held out, and counts when probit gains them. Building probit counts now
is the alternative, an engine change of roughly 25 to 30 agent-hours with its review.

What does not move. No hurdle fit that runs today changes a draw. The fit object gains only what the new inputs need, and
fitted's third argument becomes sample, as you ruled; hurdle fits are new in 1.0-0, so nothing shipped breaks.

What it costs: about 2,200 lines, nearly all R and tests, in six slices; about 45 agent-hours as estimated, 70 to 90 by the
record of plans here.

## Goal

After it lands, a hurdle.lognormal fit takes the formula interface, a data object, test, subset, weights, offset and
offset.test, with the doors reading them as on other families: extract and fitted take sample = "test", predict takes an
offset in its reserved slot and 0/1 weights, residuals and the partial dependence functions read the fit's rows, weights
and per-part offsets. Each part remains an ordinary probit or gaussian fit equal draw for draw to the standalone fit, and no
call taken today moves a draw.

## Context

- Design: docs/design/hurdle.md (the composition, section 13 authoritative); docs/plans/archive/hurdle.md (how it was built).
- R, fit: the front door's hurdle arm in [`bart`](../../R/bart.R) (the refusals of a forest() term, weights, subset, a
  variance forest, offset, offset.test, an offset() term and test, and `checkFamilyUnsupportedArgs`'s keepTrainingFits and
  keepFits refusals), [`bart2Hurdle`](../../R/bart.R), [`splitHurdleResponse`](../../R/bart.R),
  [`refuseHurdlePositiveMissingness`](../../R/bart.R), [`applyNaActionToXY`](../../R/data.R),
  [`nameOmittedRows`](../../R/data.R), [`redirectCall`](../../R/utility.R) as the components are called; the survival formula
  route in [`dbarts`](../../R/dbarts.R) (ingest once into a data object, then derive), the shape S1 follows.
- R, methods: the bartHurdle section of R/generics.R: [`hurdleParts`](../../R/generics.R),
  [`codedRowDraws`](../../R/generics.R), [`combineHurdleChannel`](../../R/generics.R), [`finishHurdle`](../../R/generics.R),
  [`hurdleLogLik`](../../R/generics.R), `extract.bartHurdle`, `fitted.bartHurdle`, `residuals.bartHurdle`,
  `predict.bartHurdle`, `print.bartHurdle`; and the shared [`trainingRowDraws`](../../R/generics.R),
  [`withTrainingOffset`](../../R/generics.R), [`refusePredictOffsetChannel`](../../R/generics.R),
  [`predictForeignReasons`](../../R/generics.R) and [`fittedForeignReasons`](../../R/generics.R) (their "sample" and
  "weights" reasons name the hurdle's absent channels), [`preparePredictRows`](../../R/data.R),
  [`padOmittedRows`](../../R/data.R), [`padPredictedRows`](../../R/data.R).
- R, offsets: [`evaluateOffsetArgument`](../../R/data.R), [`predictTermOffset`](../../R/generics.R),
  [`formulaTermOffset`](../../R/data.R), the data object's `offset.argument` attribute;
  [`refuseWithoutTestSubjects`](../../R/bart.R) for a fit with no test set.
- R, partial dependence: [`pdbart.rowSampler`](../../R/partialDependence.R), [`pdbart.frame`](../../R/partialDependence.R)
  (drops 0-weight and masked rows), [`pdbart.storedOffset`](../../R/partialDependence.R),
  [`pdbart.gridOffset`](../../R/partialDependence.R), [`pdbart.refuseNewdataOffset`](../../R/partialDependence.R),
  [`pdbart.trainingRows`](../../R/partialDependence.R), [`pdbart.buildResult`](../../R/partialDependence.R).
- Help: man/bart.Rd (the family item's hurdle paragraph, the offset item's bartHurdle sentence, the sample item, the value
  section's bartHurdle paragraph), man/pdbart.Rd where it names hurdle; inst/NEWS.Rd's hurdle bullet.
- Tests and gates: inst/tinytest/test-hurdle.R, test-hurdle-surface.R, and the hurdle arms of test-generics-errors.R,
  test-argument-surface.R, test-formula-offset-term.R, test-predict-training-rows.R, test-pdbart-averaging.R,
  test-row-names.R; benchmarks/R/hurdle-reduction.R, benchmarks/R/hurdle-exact.R (in exact-gates.yaml), the hurdle scenario
  of benchmarks/R/equivalence.R.
- TODO: hurdle feature completion; weighted-binary (docs/plans/weighted-binary.md, integer-weight probit, "scheduled after
  1.0-0 (VD 2026-08-28)").

## The rulings that stand

Every ledger entry naming hurdle, two-part, test channel, subset or offset on a composed family was read; all stand. Those
this plan builds on: dec-A15, dec-A114, dec-B07, dec-B236, dec-A19, dec-B82, dec-B14, dec-A105, dec-A112, dec-A113,
dec-A143 and dec-B379 (a missing response is dropped under the default na.action), dec-A84 (na.exclude pads on every
class), dec-A77, dec-A78, dec-A68 and dec-A147 (k, varcount and leaf.prior.sd as a list named zero and positive), dec-B207,
dec-B215, dec-B227, dec-A153, dec-B209 and dec-B377 (partial dependence: the mean response by default on a hurdle fit, the
link being "the positive part's predictor alone"), dec-B123 and dec-A102 (the residual prior on the family object), dec-B154,
dec-B184, dec-B340, dec-B341, dec-A111, dec-A240 (test weights read from the test frame), dec-B351, dec-B97, dec-B238 (a
leaf is empty only when no row reaches it), dec-B563 to dec-B565, dec-B569, dec-B570, dec-B580, dec-B13. Unaffected and
standing: dec-A25, dec-A136 and dec-B375 (keepCall), dec-A230 (the per-chain ppd layout), dec-B128 (no tombstone for
twopart), dec-B188, dec-B226. Two record the build as it stands rather than rule: dec-A225 ("multinomial and hurdle take no
weights") and dec-A228 ("hurdle has no such channel"); they describe main until this lands.

## What the two-part model defines

The hurdle's likelihood factors: sum over all rows of the zero part's Bernoulli log-likelihood of 1{y > 0}, plus sum over
the positive rows of the lognormal log density of y. The parts share no parameter, so the posterior is the product of a
probit fit on every row and a gaussian fit of log y on the positive rows, which is what the package composes.

- Test set: rows at which both parts are evaluated, nothing more. A user giving test means draws at those rows, as on every
  family.
- Subset: which rows exist. One subset for both parts (pscl's hurdle and countreg's build one model frame; lm and glm apply
  subset before na.action, run in p3-baseR.R: lm(subset = 2:10) with a missing value at row 3 gives 9 residuals, its
  na.action naming row 3).
- Weights: each row's log-likelihood times w in both parts (pscl's weights equal the expanded rows in both parts to 8e-9,
  p2-pscl.R; brms multiplies the row's whole hurdle log density by its weight and draws one observation per row). A weight
  of 0 takes a row out of both likelihoods and keeps it in both designs.
- Offset: a known shift of a part's linear predictor. The parts' predictors are on two scales, the probit latent and log y,
  so an offset belongs to a part. pscl and countreg put the offset argument in the count part only and an offset() term in
  the part whose formula carries it, a one-part formula's term in both parts (p2-pscl.R); brms puts an offset in the
  distributional parameter whose formula carries it.

## Design

### What each part's sampler receives

With n rows after subset and na.action, P the positive rows among them, and m test rows:

| | zero part (probit) | positive part (gaussian on log y) |
|---|---|---|
| x | the n rows | the rows in P |
| y | 1{y > 0} | log y |
| weights | the n weights | the weights of P |
| offset | the zero part's offset at the n rows | the positive part's offset at P |
| test | the m test rows | the n rows, then the m test rows |
| offset.test | the zero part's offset at the m rows | the positive part's offset at the n rows, then at the m rows |
| seed | the first derived seed, as today | the second, as today |

With no test set the positive part's test is the n rows, as today, so a fit taken today is unchanged. Each part is called
through bart as today, so each is the standalone fit and hurdle-reduction.R proves it draw for draw. p4-components.R ran the
premise that a part's draws do not depend on its test rows: a gaussian part fit with test = x and with test = rbind(x, xt)
gave identical training draws, sigma and test draws at the first n rows, and a probit fit with and without a test set
identical training draws.

### Ingestion (S1)

The hurdle arm ingests once, as the survival formula route does: the formula or matrix, data, test, subset, weights,
offset, offset.test and na.action are evaluated once (dec-B351) into one data object; the response checks run on its y;
then each part's data is derived from it by rows, so both parts carry the same columns, factor levels and coding. The
positive part's design must code the n rows and the test rows against the same columns as the zero part's, which the
current matrix path gets for free and a formula path must keep: a level seen only on zero rows (dec-A113 keeps it) codes on
the positive part's test rows. Each part's data object carries the formula's terms and its own offset record, so predict
and the partial dependence functions re-evaluate a part's offset on new rows as they do for any fit. A data object built by
dbartsData is taken as a formula is. A sparse matrix stays refused, by name.

The matrix route keeps today's order and draws: na.action settles the rows before the split, then the split, then the
missing-on-zero-rows refusal (dec-A113), which reads the rows as the design holds them, weight-0 rows included, as dec-B238
counts switched-off rows. The response checks (non-negative, finite, at least one zero and one positive) read the design's
rows and do not consult the weights: p12-probit-weights.R ran a gaussian and a probit fit with every weight 0, and both fit.

### na.action and the stored y

subset applies first, then na.action, then the split, as model.frame orders them (p3-baseR.R). The fit stores the
natural-scale y over the n rows and the na.action record over the rows the subset kept, so fitted and residuals pad to that
count under na.exclude (dec-A84), as lm does. A missing weight or offset value is answered as a gaussian fit given the same
input answers it, on each route. The
TODO entry's note that the real channel must reduce the stored y by the zero part's na.action is met by the order above:
no part applies an na.action of its own, so neither drops a row the hurdle's record does not.

### The test channel (S2)

- bart(..., test = ) on a hurdle fit is taken: a matrix or a frame as on other families, refused for a new missing value
  in a column complete in training (dec-A112, raised by the zero part, whose training rows are all n).
- extract(sample = "test") and fitted(sample = "test") read the zero part's test channel and the positive part's test
  columns after the first n, glued as the training rows are, every type ("ev", "ppd", "prob", "bart"); "loglik" stays
  training only, as today. fitted.bartHurdle becomes (object, type, sample, ci.level, ...), the order dec-A15 gives every
  class with a test channel; hurdle fits are new in 1.0-0, so a positional ci.level third has no shipped caller.
- A fit with no test set refuses sample = "test" in the words the other families use (`refuseWithoutTestSubjects`).
- keepTrainingFits = FALSE is taken: its refusal's reason, that the combined mean is built from the two training fits, no
  longer holds where a test set is given. The zero part drops its training channel and the positive part's test holds the m
  test rows alone; training reads are refused as on other families. keepFits = FALSE (a callback) stays refused, a callback
  having no single sampler to report.
- The fit records its test rows' names and, from a test frame's weights column (dec-A240), their weights.

### Subset and weights (S3)

- subset: an index, a logical or an expression over the data, evaluated once and applied before the split; both parts
  receive the same rows. Taken on the matrix route as bart takes it on other families (p4b-components.R ran
  bart(x, y, subset = 1:100) on a gaussian fit).
- weights: per open call 2, recommended 0 and 1 on both parts, which mean the same under the precision and frequency
  readings (dec-B564); any other value is refused by name, saying that a hurdle fit takes the weights its probit part takes.
  p12-probit-weights.R: a probit fit refuses 2 and 0.5 ("probit models do not support weights other than 0 and 1"), takes
  0/1. The positive part receives the positive rows' weights, and at weight 0 a row is out of both likelihoods and in both
  designs.
- A 0 row: fitted and residual reported, as lm reports them (p3-baseR.R); NA in the predictive draws, by the hurdle's own
  composition, which sets the row's draws to NA after drawing so the stream is unchanged (dec-B580, as dec-B472 does); NaN
  in the log-likelihood (dec-B570). p8 in p4-components.R ran today's probit at a 0 row: a finite predictive draw, which
  frequency-weights' S5 turns to NA; the hurdle reads its own record of the weights, not the zero part's, so it is NA here
  whichever lands first.
- Test weights (from a test frame): 0 and 1, a 0 giving NA predictive draws at that test row.
- predict gains weights after offset, in predict.bart's place (dec-B07), taking 0 and 1: weights = 1 draws a held-out row as
  a new observation, the door dec-B580 names, and 0 holds a row out as NA; the types that draw nothing read no weight, as on
  other families. Hurdle fits are new in 1.0-0, so combineChains moving one place breaks no shipped call. The ordinal,
  negative binomial and multinomial predicts also lack the formal dec-B580's door needs; that is frequency-weights' surface,
  noted under its interaction below.

### Offsets (S4)

Per open call 1, recommended:

- A plain offset, the argument or an offset() term, enters the positive part; the two are summed as on every family
  (dec-B154).
- offset = list(zero = , positive = ) gives each part its own; either element may be left out. The names are the ones the
  fit's k and varcount lists carry (dec-A77). Each element is a vector or an expression over the data, evaluated once
  (dec-B351); an element written in the call is recorded as an expression so predict can re-evaluate it, one handed over in
  a variable as values. A list with other names, or an element that is not numeric, is refused naming the two.
- offset.test takes the same forms and defaults to offset (bart's formal), each element evaluated on the test frame or
  reused when the rows match, refused otherwise asking for offset.test, the rule each family applies today (p10-offset-test.R:
  a formula fit's log(t) is evaluated on the test frame; a matrix fit's vector is reused for a test of n rows and refused for
  5).
- What each type carries: "prob" the zero part's offset, "bart" (the link) the positive part's, "ev" and "ppd" both
  (p4-components.R: type "bart" carries the offset on gaussian and probit fits).
- predict's reserved slot (dec-A114) takes the same forms. At newdata each part's recorded offset is re-evaluated and a given
  offset adds to its part (dec-B184); with no newdata a given offset replaces its part's at the training rows, read off the
  stored draws (dec-B340, dec-B341), by each part's `withTrainingOffset`. A part the list leaves out keeps its own. A plain
  vector the fit was given cannot be re-evaluated on newdata and is refused asking for offset, per part.
- pdbart and pd2bart carry each part's stored offset to predict as the list, and refuse a newdata a part's offset cannot be
  evaluated on, naming the part.

### Every door

| door | today | after |
|---|---|---|
| bart(formula / data object) | refused | taken (S1) |
| bart(test, subset, weights, offset, offset.test) | refused | taken (S2 to S4); weights 0/1 |
| bart(keepTrainingFits = FALSE) | refused | taken (S2) |
| bart(keepFits = FALSE, callback, samplerOnly, warm.start, variance) | refused | refused, unchanged |
| dbarts(), dbartsSpec, xbart | refused | refused, unchanged (dec-B82, dec-B14) |
| extract(sample = "test") | refused | taken (S2) |
| fitted | (object, type, ci.level) | (object, type, sample, ci.level) (S2, dec-A15) |
| residuals | training rows | training rows, the 0 rows reported (S3) |
| predict(offset) | refused | taken (S4) |
| predict(weights) | no formal | 0 and 1 (S3) |
| pdbart, pd2bart | the rows of the fit, no offset | the fit's rows less its 0 rows, each part's offset (S3, S4) |

### What changes on the fit object

- `$zero` and `$positive` carry their weights, offsets and test channels as any bart fit does; `$positive`'s test channel is
  n + m rows wide.
- The hurdle fit records its weights and test weights where given, its test rows' names, and keeps y, na.action and the
  training rows' names as today. Nothing a fit taken today carries is removed or renamed.

### Calls made

1. The formula interface and a data object are in scope: subset as an expression, weights and offsets as columns, test
   frames and an offset recorded for predict all need it, every other family takes it (p11-formula-families.R), and
   dec-B97 finished the survival families' interface on the same ground.
2. keepTrainingFits = FALSE is taken once a test channel exists; keepFits = FALSE stays refused.
3. predict gains weights after offset, taking 0 and 1, the door dec-B580 names.
4. The response and missingness checks read the design's rows, weights not consulted.
5. The positive part's test is the n training rows stacked over the test rows, one channel, rather than a second fit or a
   replay of trees, so extract answers without saved trees.

## Open calls

### 1. Which part a plain offset enters

Background. The zero part's predictor is a probit latent and the positive part's is log y: an offset of log exposure on the
positive part makes a positive amount proportional to exposure, while on the probit latent the same number has no such
meaning. The package refused one leaf-prior spread on a hurdle fit for that reason (dec-A105), and reports the positive
part as a hurdle fit's link (dec-B207). pscl's and countreg's offset argument enters the count part only; their offset()
term in a one-part formula enters both parts and, in a two-part formula, the part it is written in (p2-pscl.R). brms enters
an offset in the parameter whose formula carries it, the size for the main formula. On this package a term and the argument
behave alike (dec-B154, dec-B184), so pscl's split between them is not open.

- (a) A plain offset, argument or term, enters the positive part; offset = list(zero = , positive = ) states each part's.
  Recommended. About 220 lines of R and 350 of tests.
- (b) A plain offset enters both parts, as pscl's one-part term does; the list still states each. One number then reads in
  two units, which dec-A105 refused for the spread. About the same cost.
- (c) A plain offset is refused on a hurdle fit; only the named list is taken. Unambiguous, but a user's log exposure, which
  pscl, countreg and brms all put in the size part, is refused. About 20 lines fewer than (a).

Evidence that would change it: a user population whose offsets belong in the zero part, such as a known baseline
probability, for which (a) still works through the list.

### 2. What weights a hurdle fit takes

Background. A weight in a two-part model is read by pscl and brms as a count of identical rows in both parts, which
frequency weights are (dec-B563). Under that reading the zero part needs probit counts: a probit row of weight w is w latent
draws, exact for whole w only. Probit takes 0 and 1 alone today (p12-probit-weights.R), integer-weight probit is parked in
docs/plans/weighted-binary.md as scheduled after 1.0-0 (VD 2026-08-28), and dec-B563 kept probit out of frequency weights
without scheduling it. A precision reading, noise variance over w on log y, has no meaning on the probit part (dec-B564
refuses it there), so a weight read as a precision on one part and ignored on the other is not offered. Survey or
inverse-probability weights are neither (dec-B568's caution).

- (a) 0 and 1 on both parts now, the membership a probit takes, read the same under both readings; whole numbers when
  integer-weight probit is built, on probit and the hurdle together. Recommended: the hurdle's zero part stays a fit a plain
  probit call reproduces, so fitting the two parts separately remains the same posterior, as the help promises. About 100
  lines of R and 250 of tests (S3).
- (b) Build integer-weight probit now, for probit and the hurdle together, and take whole-number weights on the hurdle. An
  engine slice in the draws-changing tier (the latent update draws w truncated normals per row; tests/cpp, an exact gate, an
  opus review with mutants), moving the 2026-08-28 schedule and adding probit to dec-B563's families; roughly 25 to 30
  agent-hours with review over (a), after which frequency-weights' resolution table gains probit.

Evidence that would change it: users of semicontinuous responses with aggregated rows, whose counts (a) cannot take.

## Constraints

- No engine, bridge, header, state-format or `*.in` file changes; `git diff --stat` on the branch names none.
- Every hurdle call taken today keeps its draws: the equivalence hurdle scenario bitwise, hurdle-reduction.R bitwise,
  hurdle-exact.R at its tolerance, test-hurdle.R unchanged.
- Each part stays the standalone probit or gaussian fit at its derived seed, draw for draw, under every new input.
- Each argument evaluated once (dec-B351); both parts see the same rows.
- Out of scope: whole-number or real weights (open call 2, (a)); a variance forest on the positive part; keepFits = FALSE and
  callbacks; dbarts(), dbartsSpec, xbart and the sampler doors (dec-B14); a sparse design; a two-part formula giving each part
  its own predictors; an engine hurdle (dec-B82).
- Base R calls within DESCRIPTION's R floor; no new dependency.

## Slices

Each slice is a commit series on hurdle-completion; the feature merges once after S6.

### S1: ingestion and the formula interface (~200 R, ~30 help, ~250 tinytest)

The hurdle arm ingests once into a data object and derives the parts; the formula, a data object and a frame taken; the
matrix route's order kept. Pins: a formula fit equals the matrix fit on the same coded design draw for draw; factor columns,
a level seen only on zero rows (the positive part's test rows code it); row names; na.action records; the matrix fit of
test-hurdle.R unchanged.

### S2: the test channel (~140 R, ~40 help, ~300 tinytest)

The stacked positive test, extract and fitted at sample = "test", fitted's formals, keepTrainingFits = FALSE, the test rows'
names and weights. Pins: extract(sample = "test") at types "prob", "bart" and "ev" equals predict(newdata = test) on a fit
kept with its trees, to rounding (the stored draws against the replayed trees), these types drawing nothing; the training channel unchanged by a test set; the "ppd" test
draws' NA at a 0 test weight; refusal words without a test set.

### S3: subset and weights (~100 R, ~25 help, ~250 tinytest)

subset on both routes; 0/1 weights; predict's weights formal; the 0 row's NA draws and NaN log-likelihood; partial
dependence leaving 0 rows out. Pins: each part equals the standalone fit with the same subset and weights; a subset's
na.action record against lm's on the same rows; residuals at a 0 row.

### S4: offsets (~220 R, ~40 help, ~350 tinytest)

The plain and named forms, terms, offset.test, predict's reserved slot at newdata and with no newdata, partial dependence.
Pins: "prob" moves with a zero-part offset only and "bart" with a positive-part one only; predict at the training rows as
newdata equals fitted where offsets are expressions; predict with no newdata and the fit's own offset equals fitted;
refusals naming the part.

### S5: gates (~200 in benchmarks/R)

hurdle-reduction.R gains arms for a formula fit, a test set, a subset, 0/1 weights and each part's offset, each part
bitwise against its standalone fit built by hand. hurdle-exact.R gains an arm with a 0 row in each cell and one with an
offset in each part, the expected posterior derived from the response, weights and offsets alone (the split re-derived, as
its header requires). The equivalence hurdle scenario is not re-recorded.

### S6: records (~40)

Ledger entries for the calls made and the two rulings; this plan's Status and Landing note; TODO; the plans index; NEWS's
hurdle bullet (hurdle is new in 1.0-0, so one bullet against 0.9-34, which had no hurdle family, run in scratch/libs/cran).

## Tests and gates

Per slice, against its own library (`R CMD INSTALL -l <lib> .`, `buildInfo()$mode` checked): the touched test files singly;
the four snapshot files on the reference build; `lintr::lint_package()`, `air format --check .`, check-rc-codoc,
check-win-drift, check-doc-freshness, each on its own exit status. At the feature's end: the full tinytest suite in one
process at its default ncpu; `R CMD check --as-cran` from a clean tarball; the exact gates of exact-gates.yaml with `quick`
(the fit object changes); hurdle-reduction.R; the equivalence compare counted per scenario, the hurdle scenario among the
identical; CI on the branch tip. No tests/cpp or sanitizer run is owed: no compiled code moves.

## Mutants

A reviewer applies each to a copy, reinstalls, and runs the named killer; each must fail it.

1. S1: the positive part coded against its own rows' factor levels - the zero-rows-only level pin.
2. S1: na.action applied after the split - the matrix fit's draws in test-hurdle.R and the na.action record pin.
3. S2: the test columns read from column n, not n + 1 - the extract against predict pin.
4. S2: the positive part's test without the n training rows when a test set is given - the training channel pin.
5. S3: the positive part given the first |P| weights rather than P's - hurdle-reduction.R's weights arm and hurdle-exact.R's
   0-row arm.
6. S3: the 0 row's predictive draw left drawn - the NA pin, training and test.
7. S3: subset applied after na.action - the record pin against lm.
8. S4: a plain offset sent to both parts - the "prob" pin.
9. S4: the zero element sent to the positive part - hurdle-exact.R's offset arm.
10. S4: with no newdata, a given offset added rather than replacing - the predict-equals-fitted pin.
11. S4: partial dependence passing the zero part's offset alone - a pdbart "prob" pin with a positive-part offset (no
    change) and a zero-part one (a change).
12. Any slice: the matrix route's seed derivation or test stacking changed - the equivalence hurdle scenario.

## RNG class

Neutral. No call taken today moves a draw: today's calls are the matrix route with none of the new inputs, and the positive
part's test with no test set stays the n rows. fitted's third positional argument changes meaning on hurdle fits only, a
1.0-0 family. New calls draw what the standalone parts draw.

## Consumers

stan4bart (27ef163), bartCause (1e4b882), treatSens (268b2db) and bairrtt (6575ba0) do not name hurdle anywhere in their R
or src trees (grep, 2026-11-03); none changes.

## Interaction with frequency-weights, and order

Files both touch: R/bart.R (frequency-weights adds weights.type to bart's formals and forwarding; this plan rewrites the
hurdle arm and bart2Hurdle), R/generics.R (frequency-weights edits sampleFromPPD, pointwiseLogLikelihood and the bart
methods; this plan the bartHurdle section and the foreign-argument reason tables), man/bart.Rd (the weights item and a new
weights.type item there; the hurdle paragraph and the value section here), inst/NEWS.Rd, and test-zero-weights.R if either
adds a hurdle arm. The functions are disjoint; the conflicts are textual.

Order recommended: this feature first. It is R only and need not wait on coefficient-prior-live, which frequency-weights
follows. frequency-weights then amends two things in its S5: its resolution table's hurdle row becomes probit's (NULL
resolving to frequency membership, "precision" refused, "frequency" 0 and 1 only), and the hurdle arm routes weights.type to
both parts, refusing "precision" in the hurdle's words before the zero part refuses it in probit's; about 20 lines and a
test. The hurdle reads its own weights record, so frequency-weights' S4, which keeps probit's 0/1 weights on the data object
instead of the active-row mask, does not move it, and its own NA composition at 0 rows does not wait on S5's. If
frequency-weights lands first instead, S3 here writes weights in weights.type's terms and rebases over its R/bart.R and
man/bart.Rd edits. Either way the feature that lands second states the hurdle's weights.type row.

dec-B580 names predict(weights = 1) as the door to a 0 row's draw; the ordinal, negative binomial and multinomial predicts
have no weights formal on main or in the frequency-weights plan, which is that plan's to settle.

## Not verified

- That a formula route deriving both parts from one data object reproduces the matrix route's draws for a numeric design;
  S1's first pin shows it or names what differs.
- The exact cost of the n + m wide positive test channel at large n; it is n wide today, and m more columns are what any
  bart fit's test set costs.
- How a missing weight or offset value is answered on the matrix route, where the hurdle applies na.action itself; S3 and S4
  match a gaussian fit given the same input and pin it.
- The cost of integer-weight probit in open call 2 (b) is scaled from frequency-weights' Polya-Gamma slice and the
  weighted-binary memo, not planned.
- hurdle-exact.R's configuration was read for its header only; S5 confirms a 0 row and an offset fit within its two cells.
