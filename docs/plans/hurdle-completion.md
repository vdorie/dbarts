# hurdle-completion: a test set, weights, subset, offsets and the formula interface on hurdle.lognormal fits

Status: PLANNED 2026-11-03 on main (170635476); amended the same day after its blind critique (see
[The critique and what was amended](#the-critique-and-what-was-amended)). Not built. Both calls are ruled (dec-B587, dec-B588).
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
every family; "0 weight means out of sample"), dec-A89 (a missing weight or offset at fit: na.omit and na.exclude drop the
row), dec-B470 (predict with no newdata takes weights, which replace the fit's own), dec-B482 (with no newdata a missing
weight or offset follows na.action), dec-B185 (a recorded expression freezes its non-data names), dec-B166 (a multinomial
offset is a rows x K matrix). The runs behind it are the probes under scratch/hrd-plan/probes, on a shipped build of main
8fce28568 (scratch/libs/hrd-plan), 0.9-34 (scratch/libs/cran) and pscl 1.5.9 and countreg 0.3.0
(scratch/hrd-plan/landscape-lib); the critique's under scratch/hrd-crit/probes and the amendment's under
scratch/hrd-plan/probes-amend, on a shipped build of main 170635476 (scratch/libs/hrd-crit). brms 2.23.0 was read from its
installed source and glmmTMB from its help page (man/glmmTMB.Rd on the package's GitHub master, read 2026-11-03).

agent: S1 to S4 one sonnet implementer each (effort medium); S5 sonnet. One sonnet review for S1 and S2 together and one for
S3 and S4, each reviewer checking its fix round and running the before/after gate (S1) itself, not trusting the report; S5
takes the R-surface tier's review as a gate-only slice ([Process by risk](README.md#process-by-risk)). Tier "R surface, no
draw change": no engine, bridge, header or state-format file moves, and no call taken today changes a draw.
rng: neutral. Every hurdle call taken today (a matrix or frame x with none of test, weights, subset, offset) keeps its draws
and every method's output bit for bit, shown by the before/after gate (S1) against a library built from this plan's base;
the equivalence harness's hurdle scenario, hurdle-reduction.R, hurdle-exact.R and test-hurdle.R's pins do not move either.
New calls take new draws, each part held draw for draw to the standalone probit or gaussian fit it is.
window: before 1.0-0. The basis is the maintainer's note of 2026-09-26 in the TODO entry, "part of feature completion, not
to be forgotten"; the entry also says the schedule is not yet set, and dec-A15 records the note, not a schedule. The
work changes public surface, fitted.bartHurdle's third argument among it, which a release would then owe compatibility.
One branch, hurdle-completion, merged once (dec-B494). R only, so it need not wait on coefficient-prior-live or the engine
queue; order against frequency-weights below.
budget: about 2,600 lines: R ~710, help ~145, NEWS ~5, tinytest ~1,350, gate scripts ~350, records ~50. Agent time,
estimated and not measured: S1 ~7.5, S2 ~7, S3 ~6, S4 ~9.5, S5 ~5.5, S6 ~1, three reviews with fix rounds ~12, gates and
check ~3.5: ~52. Plans here run 1.5 to 2 times low, so 80 to 105.

## Summary for the maintainer

A hurdle fit today takes only predictors as a matrix or data frame, and a response. It refuses a formula, a test set,
weights, a subset and an offset, and its fitted method has no way to ask for test rows. It is the only family that refuses
a formula. This plan finishes it so that a hurdle fit is made the way any other fit is made.

What a user gets. bart(y ~ ., data, family = "hurdle.lognormal") works, with subset, weights, an offset and a test set
written as on any other family. The test set gives stored draws at those rows, so extract and fitted take sample = "test"
and a fit kept without its trees still answers for them. Both parts always see the same rows: a subset or a missing value
removes a row from both, as one model frame does in lm and pscl's hurdle. predict takes an offset in the place held for it
and, as a gaussian or probit fit's predict does, a weight of 1 to draw a held-out row as a new observation.

Two choices were the maintainer's, and both are ruled as recommended (see [Ruled](#ruled)).

First, which part an offset enters. The zero part is a probit and the positive part models log y, so one number cannot
mean the same thing in both, the reason a single leaf-prior spread is refused on this family. pscl and countreg put their
offset argument in the count part, and brms puts an offset in the part whose formula carries it, the size part for the
main formula; glmmTMB's help says its offset argument is "for conditional model (only)". The package already reports the
positive part as a hurdle fit's link scale. I recommend that a plain offset enters the positive part, and that a list
naming zero and positive gives each part its own, the names extract already gives this family's k and variable counts. No
argument takes a per-part value today, so the list is a new shape; the one the package has is the multinomial's offset, a
matrix with one named column per category, the alternative here. One difference from pscl the help must say plainly: in
pscl, y ~ x + offset(log(t)) reaches both parts, and here it reaches the positive part only.

Second, what weights it takes. pscl and brms read a hurdle's weights as counts of identical rows in both parts. The zero
part is a probit fit, and a probit fit takes only 0 and 1, which say whether a row is in the data; counts on probit were
scheduled after 1.0-0 in August and left unscheduled when you ruled frequency weights. I recommend 0 and 1 now, read as you
ruled for every family on 2026-11-03 (0 holds a row out, and predict with a weight of 1 draws it as a new observation), and
counts when probit gains them; a user with counts expands the rows meanwhile. That slice, subset with it, is about 6
agent-hours. Building probit counts now is the alternative, an engine change of roughly 25 to 30 agent-hours more with its
review. Ruled: option (a) of each.

One call made for you to know of: a missing weight or offset on a hurdle fit's matrix interface is dropped under na.omit
and na.exclude, as your 2026-09-27 rule says fitting does; gaussian and probit fits there still refuse it whatever
na.action says, a gap filed for later.

What does not move. No hurdle fit that runs today changes a draw, and a script run on the same machine against the
package before and after shows every such call's draws and outputs identical. The fit object gains only what the new inputs
need, and fitted's third argument becomes sample, as you ruled; hurdle fits are new in 1.0-0, so nothing shipped breaks.

What it costs: about 2,600 lines, nearly all R and tests, in six slices; about 52 agent-hours as estimated, 80 to 105 by
the record of plans here.

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
  route in [`dbarts`](../../R/dbarts.R) (ingest once into a data object, then derive), the shape S1 follows;
  [`dbartsData`](../../R/data.R)'s warning that arguments beside a data object are ignored.
- R, methods: the bartHurdle section of R/generics.R: [`hurdleParts`](../../R/generics.R),
  [`codedRowDraws`](../../R/generics.R), [`combineHurdleChannel`](../../R/generics.R), [`finishHurdle`](../../R/generics.R),
  [`hurdleLogLik`](../../R/generics.R), `extract.bartHurdle`, `fitted.bartHurdle`, `residuals.bartHurdle`,
  `predict.bartHurdle`, `print.bartHurdle`, `summary.bartHurdle`, and [`plot.bartHurdle`](../../R/plot.R) (its log-scale
  panel selects the positive part's test columns with a length-n logical); and the shared
  [`trainingRowDraws`](../../R/generics.R), [`withTrainingOffset`](../../R/generics.R) (its ordinary branch shifts the
  training channel, which on the positive part is the positive rows only), [`shiftObservations`](../../R/generics.R),
  [`trainingRowsAnswered`](../../R/generics.R), [`fillMissingChannel`](../../R/generics.R),
  [`hazardSurvivalProbabilities`](../../R/bart.R) (the precedent for shifting a stored test channel against offset.test with
  no newdata), [`refusePredictOffsetChannel`](../../R/generics.R),
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
  section's bartHurdle paragraph), man/dbarts.Rd's hurdle item (its sentence that weights, subset, offset and test "are all
  unsupported this arc"), man/pdbart.Rd where it names hurdle; inst/NEWS.Rd's hurdle bullet.
- Tests and gates: inst/tinytest/test-hurdle.R, test-hurdle-surface.R, and the hurdle arms of test-generics-errors.R,
  test-argument-surface.R, test-formula-offset-term.R, test-predict-training-rows.R, test-pdbart-averaging.R,
  test-row-names.R; benchmarks/R/hurdle-reduction.R (whose standalone fits today take their split from the package's own
  splitter), benchmarks/R/hurdle-exact.R (in exact-gates.yaml), the hurdle scenario of benchmarks/R/equivalence.R; a new
  before/after script, benchmarks/R/hurdle-before-after.R (S1).
- TODO: hurdle feature completion; weighted-binary (docs/plans/weighted-binary.md, integer-weight probit, "scheduled after
  1.0-0 (VD 2026-08-28)"); fit-na-action-weights-offset (the matrix route's missing weight or offset, filed with this
  amendment).

## The rulings that stand

Every ledger entry naming hurdle, two-part, test channel, subset or offset on a composed family was read; all stand. Those
this plan builds on: dec-A15, dec-A114, dec-B07, dec-B236, dec-A19, dec-B82, dec-B14, dec-A105, dec-A112, dec-A113, dec-A89,
dec-B470, dec-B482, dec-B185, dec-B166,
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
  the part whose formula carries it, a one-part formula's term in both parts (p2-pscl.R, c05-pscl.R); brms puts an offset
  in the distributional parameter whose formula carries it; glmmTMB's help gives its offset argument to the conditional
  (non-zero) model only and, when the zero formula copies the main one (~ .), drops its offset() term ("the offset term
  will automatically be dropped").

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
and the partial dependence functions re-evaluate a part's offset on new rows as they do for any fit. Under the ruling on which part a plain offset enters,
(a) the zero part's terms drop the formula's offset() term, which is the positive part's, so neither predict(fit$zero,
newdata) through [`predictTermOffset`](../../R/generics.R) nor [`pdbart.storedOffset`](../../R/partialDependence.R) and
[`pdbart.gridOffset`](../../R/partialDependence.R), which read the zero part's sampler, add it to the probit; the zero part's
offset record holds the zero element alone. A data object built by dbartsData is taken as a formula is, and as on every
family it is the whole data: an argument beside it, offset, weights, subset or test, is ignored with the package's warning
("if data supplied as dbartsData, remaining arguments are ignored"; a01-amend.R, on a gaussian fit). Its one offset slot is
a plain offset and enters the part call 1 sends a plain offset to; bart(dataObject, offset = list(zero = , positive = ))
therefore warns that offset is ignored, and a per-part offset is written through the formula. A sparse matrix stays
refused, by name.

The matrix route keeps today's order and draws: na.action settles the rows before the split, then the split, then the
missing-on-zero-rows refusal (dec-A113), which reads the rows as the design holds them, weight-0 rows included, as dec-B238
counts switched-off rows. The response checks (non-negative, finite, at least one zero and one positive) read the design's
rows and do not consult the weights: p12-probit-weights.R ran a gaussian and a probit fit with every weight 0, and both fit.

### na.action and the stored y

subset applies first, then na.action, then the split, as model.frame orders them (p3-baseR.R). The fit stores the
natural-scale y over the n rows and the na.action record over the rows the subset kept, so fitted and residuals pad to that
count under na.exclude (dec-A84), as lm does. A missing weight or offset value (either part's) marks its row incomplete, and
na.action answers it on both routes as dec-A89 states fitting does: na.omit and na.exclude drop the row from both parts,
the default and na.pass refuse it. A gaussian fit's formula route does this through its model frame (c04b-formula-na.R),
and the hurdle's formula route inherits it. On the matrix route a gaussian or probit fit refuses under every na.action
("'weights' cannot be NA", "'offset' contains missing values"; c04b-formula-na.R; 0.9-34 refused too,
a03-cran-matrix-na.R); the hurdle's matrix route applies na.action itself, so it folds the weight and offset columns into
the rows its na.action step reads and follows dec-A89 directly. That gap on the other families is filed in TODO
(fit-na-action-weights-offset), not closed here. The TODO entry's note that the real channel must reduce the stored y by
the zero part's na.action is met by the order above: no part applies an na.action of its own, so neither drops a row the
hurdle's record does not.

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
  test rows alone, so test reads start at column 1 there and at n + 1 otherwise. Its test draws equal the default fit's at
  one seed bit for bit (a02-test-tail.R: a gaussian part's draws at m test rows are identical whether its test holds them
  alone or after the n training rows, at one and two threads). Training reads are refused in the words a gaussian fit kept
  so uses (c07-weights.R, a01-amend.R): fitted and residuals ("cannot extract train sample predictions; bart must be called
  with 'keepTrainingFits' == TRUE ..."), extract(type = "loglik") and predict with no newdata alike, and plot ("plot requires
  bart to be called with 'keepTrainingFits' == TRUE ..."); summary, print, predict at newdata and the partial dependence
  functions answer, as on a gaussian fit. keepFits = FALSE (a callback) stays refused, a callback having no single sampler
  to report.
- plot.bartHurdle's log-scale panel reads the first n columns of the positive part's test channel; today it selects with a
  length-n logical, which over n + m columns recycles and fails with "'x' and 'y' lengths differ" whenever a positive row is
  among the first m rows (c10-plot-stacked.R).
- The fit records its test rows' names and, from a test frame's weights column (dec-A240), their weights.

### Subset and weights (S3)

- subset: an index, a logical or an expression over the data, evaluated once and applied before the split; both parts
  receive the same rows. Taken on the matrix route as bart takes it on other families (p4b-components.R ran
  bart(x, y, subset = 1:100) on a gaussian fit).
- weights: per the weights ruling (dec-B588), 0 and 1 on both parts, which mean the same under the precision and frequency
  readings (dec-B564); any other value is refused by name, saying that a hurdle fit's weights mark rows in (1) and out (0)
  and that whole-number counts are fit by repeating each row its weight's number of times, which the probit part reads
  exactly. p12-probit-weights.R: a probit fit refuses 2 and 0.5 ("probit models do not support weights other than 0 and
  1"), takes 0/1. The positive part receives the positive rows' weights, and at weight 0 a row is out of both likelihoods
  and in both designs.
- A 0 row: fitted and residual reported, as lm reports them (p3-baseR.R); NA in the predictive draws, by the hurdle's own
  composition, which sets the row's draws to NA after drawing so the stream is unchanged (dec-B580, as dec-B472 does); NaN
  in the log-likelihood (dec-B570). p8 in p4-components.R ran today's probit at a 0 row: a finite predictive draw, which
  frequency-weights' S5 turns to NA; the hurdle reads its own record of the weights, not the zero part's, so it is NA here
  whichever lands first.
- Test weights (from a test frame): 0 and 1, a 0 giving NA predictive draws at that test row.
- predict gains weights after offset, in predict.bart's place (dec-B07), taking 0 and 1. The rule rests on dec-B580's doors,
  not on the probit part's predict, which takes positive whole trial counts and refuses 0 (c07-weights.R): weights = 1
  draws a held-out row as a new observation and 0 holds a row out as NA; the types that draw nothing read no weight, as on
  other families. With no newdata, given weights replace the hurdle's own for the training rows' predictive draws, one value
  or one per training row (the hurdle's n, not the positive part's), worked from the stored draws, as dec-B470 rules for
  every family; a missing value there follows na.action (dec-B482). Hurdle fits are new in 1.0-0, so combineChains moving
  one place breaks no shipped call. The ordinal, negative binomial and multinomial predicts also lack the formal dec-B580's
  door needs; that is frequency-weights' surface, noted under its interaction below.

### Offsets (S4)

Per the offset ruling (dec-B587):

- A plain offset, the argument or an offset() term, enters the positive part; the two are summed as on every family
  (dec-B154).
- offset = list(zero = , positive = ) gives each part its own; either element may be left out. The names are those
  extract gives the fit's k and varcount (dec-A77); no hurdle argument takes a per-part value today (k = list(zero = ,
  positive = ) is refused, c02-ingest.R), so the list is a new input shape. Each element is a vector or an expression over
  the data, evaluated once (dec-B351); an element written in the call is recorded as its own expression, its non-data names
  frozen and no environment kept (dec-B185), so predict can re-evaluate it, one handed over in a variable as values. A list
  with other names, or an element that is not numeric, is refused naming the two.
- The help says plainly that y ~ x + offset(log(t)) reaches the positive part only here, where pscl's hurdle puts that term
  in both parts, and that a zero-part offset is written offset = list(zero = ).
- offset.test takes the same forms and defaults to offset (bart's formal), each element evaluated on the test frame or
  reused when the rows match, refused otherwise asking for offset.test, the rule each family applies today (p10-offset-test.R:
  a formula fit's log(t) is evaluated on the test frame; a matrix fit's vector is reused for a test of n rows and refused for
  5).
- What each type carries: "prob" the zero part's offset, "bart" (the link) the positive part's, "ev" and "ppd" both
  (p4-components.R: type "bart" carries the offset on gaussian and probit fits).
- predict's reserved slot (dec-A114) takes the same forms. At newdata each part's recorded offset is re-evaluated and a given
  offset adds to its part (dec-B184). With no newdata a given offset replaces its part's at the training rows, read off the
  stored draws (dec-B340, dec-B341), one value or one per training row (the hurdle's n). The zero part's goes through
  `withTrainingOffset`, its training channel being the n rows. The positive part's cannot: the hurdle reads that part's n
  training-row values from its test channel, while `withTrainingOffset` shifts its training channel, the positive rows
  only, and refuses an offset of length n (a01-amend.R: "'offset' must be one number or one per training row (37)"). The
  positive part's replacement shifts the first n columns of its test channel by the new offset less offset.test[1:n], which
  the part stores, as `hazardSurvivalProbabilities` does for stored test subjects; a01-amend.R ran that shift against
  predict at the n rows as newdata with the same offset: equal to 4e-16. A missing value in either part's replacement
  follows na.action as dec-B482 rules (the default refusing, na.pass an NA row, na.omit dropping it), the row being one row
  of both parts. A part the list leaves out keeps its own. A plain vector the fit was given cannot be re-evaluated on
  newdata and is refused asking for offset, per part.
- pdbart and pd2bart carry each part's stored offset to predict as the list, and refuse a newdata a part's offset cannot be
  evaluated on, naming the part.

### Every door

| door | today | after |
|---|---|---|
| bart(formula / data object) | refused | taken (S1) |
| bart(frame x) | taken | taken, unchanged (S1, before/after gate) |
| bart(data object, offset / weights / subset / test beside it) | refused | ignored with dbartsData's warning, as on every family (S1) |
| bart(test, subset, weights, offset, offset.test) | refused | taken (S2 to S4); weights 0/1 |
| a missing weight or offset at fit | refused | na.action answers it, dec-A89 (S3, S4) |
| bart(keepTrainingFits = FALSE) | refused | taken (S2); training reads, residuals, loglik and plot refused in a gaussian fit's words |
| bart(keepFits = FALSE, callback, samplerOnly, warm.start, variance) | refused | refused, unchanged |
| dbarts(), dbartsSpec, xbart | refused | refused, unchanged (dec-B82, dec-B14) |
| extract(sample = "test") | refused | taken (S2) |
| fitted | (object, type, ci.level) | (object, type, sample, ci.level) (S2, dec-A15) |
| residuals | training rows | training rows, the 0 rows reported (S3) |
| plot | the n columns | the first n columns of n + m (S2) |
| predict(offset) | refused | taken; with no newdata the positive part's test channel shifted (S4) |
| predict(weights) | no formal | 0 and 1; with no newdata they replace the fit's own, dec-B470 (S3) |
| pdbart, pd2bart | the rows of the fit, no offset | the fit's rows less its 0 rows, each part's offset (S3, S4) |

### What changes on the fit object

- `$zero` and `$positive` carry their weights, offsets and test channels as any bart fit does; `$positive`'s test channel is
  n + m rows wide.
- The hurdle fit records its weights and test weights where given, its test rows' names, and keeps y, na.action and the
  training rows' names as today. Nothing a fit taken today carries is removed or renamed.
- The positive part keeps its own training channel (the positive rows), though no hurdle method reads it. Building that
  part always with keepTrainingFits = FALSE would save |P| of the 2n + |P| + 2m values a draw stores, moving no draw; it is
  not taken (call made 8).

### Calls made

1. The formula interface and a data object are in scope: subset as an expression, weights and offsets as columns, test
   frames and an offset recorded for predict all need it, every other family takes it (p11-formula-families.R), and
   dec-B97 finished the survival families' interface on the same ground.
2. keepTrainingFits = FALSE is taken once a test channel exists; keepFits = FALSE stays refused.
3. predict gains weights after offset, taking 0 and 1, the door dec-B580 names.
4. The response and missingness checks read the design's rows, weights not consulted.
5. The positive part's test is the n training rows stacked over the test rows, one channel, rather than a second fit or a
   replay of trees, so extract answers without saved trees.
6. A data object is the whole data on a hurdle fit as on every family: offset, weights, subset or test beside it are
   ignored with dbartsData's warning, and its one offset slot is a plain offset. The alternative, taking an offset list
   beside a data object on this family alone, is a one-family exception to the data object's rule; a per-part offset has the
   formula.
7. A missing weight or offset on the hurdle's matrix route follows dec-A89 (na.omit and na.exclude drop the row, the default
   and na.pass refuse it), where gaussian and probit fits on that route refuse under every na.action. The alternative was
   matching those fits, as the plan first had it, which keeps a gap against dec-A89's stated rule. To be told to the
   maintainer with the landing; the other families' gap is a TODO entry.
8. The positive part keeps its training channel. Its saving, |P| values per draw, is at most a third of what the hurdle
   stores, and keepTrainingFits = FALSE already gives the larger one, while without it fit$positive would stop being the
   ordinary gaussian fit the help says it is (its fitted, residuals, loglik and plot refusing) and hurdle-reduction.R could
   not compare its training draws with the standalone fit's.

## Ruled

### 1. Which part a plain offset enters (ruled: (a), dec-B587)

Background. The zero part's predictor is a probit latent and the positive part's is log y: an offset of log exposure on the
positive part makes a positive amount proportional to exposure, while on the probit latent the same number has no such
meaning. The package refused one leaf-prior spread on a hurdle fit for that reason (dec-A105), and reports the positive
part as a hurdle fit's link (dec-B207). pscl's and countreg's offset argument enters the count part only; their offset()
term in a one-part formula enters both parts and, in a two-part formula, the part it is written in (p2-pscl.R,
c05-pscl.R). glmmTMB's help gives its offset argument to the conditional model only and drops an offset() term from a zero
formula written as the main formula's right-hand side. brms enters an offset in the parameter whose formula carries it,
the size for the main formula. On this package a term and the argument behave alike (dec-B154, dec-B184), so pscl's split
between them is not open; under (a) a pscl user's y ~ x + offset(log(t)) reaches both parts in pscl and only the positive
part here, which the help says plainly.

A per-part offset needs a shape, and no hurdle argument takes a per-part value today (k = list(zero = , positive = ) is
refused, c02-ingest.R); the zero and positive names are those extract gives (dec-A77). The package's one per-component
offset is the multinomial's, an n x K matrix whose columns are matched to the categories by name (dec-A111; dec-B166 on the
flat C API), refusing a list or a vector (c02-ingest.R). Here that would read offset = cbind(zero = , positive = ). The list
is recommended over it: either part can be left out, where a matrix needs a column of zeros, and each element is recorded
as its own expression for predict (dec-B185), where cbind(zero = 0, positive = log(t)) is one expression whose parts
predict must split. The matrix's cost is the reverse: one shape for every per-component offset in the package.

Ruled (maintainer, 2026-11-03, dec-B587), option (a): "Ok, use your recommendation for this item." A plain offset, argument
or offset() term, enters the positive part; offset = list(zero = , positive = ) states each part's. Before the ruling no
ruling settled it: dec-A114 says only that hurdle offsets are coming, and dec-A105 is about a leaf-prior spread.

- (a) A plain offset, argument or term, enters the positive part; offset = list(zero = , positive = ) states each part's.
  Recommended. About 9.5 agent-hours (S4).
- (b) A plain offset enters both parts, as pscl's one-part term does; the list still states each. One number then reads in
  two units, which dec-A105 refused for the spread. About the same.
- (c) A plain offset is refused on a hurdle fit, a data object's offset slot with it; only the named list is taken.
  Unambiguous, but a user's log exposure, which pscl, countreg, glmmTMB and brms all put in the size part, is refused.
  About half an hour less than (a).

Evidence that would change it: a user population whose offsets belong in the zero part, such as a known baseline
probability, for which (a) still works through the list.

### 2. What weights a hurdle fit takes (ruled: (a), dec-B588)

Background. A weight in a two-part model is read by pscl and brms as a count of identical rows in both parts, which
frequency weights are (dec-B563). Under that reading the zero part needs probit counts: a probit row of weight w is w latent
draws, exact for whole w only. Probit takes 0 and 1 alone today (p12-probit-weights.R), integer-weight probit is parked in
docs/plans/weighted-binary.md as scheduled after 1.0-0 (VD 2026-08-28), and dec-B563 kept probit out of frequency weights
without scheduling it. A precision reading, noise variance over w on log y, has no meaning on the probit part (dec-B564
refuses it there), so a weight read as a precision on one part and ignored on the other is not offered. Survey or
inverse-probability weights are neither (dec-B568's caution).

Ruled (maintainer, 2026-11-03, dec-B588), option (a): "Let's start with your recommendation and think about extensions
in the future." 0 and 1 on both parts now; other values refused, naming the remedy of repeating each row its count's number
of times. Before the ruling nothing settled it: dec-B563 left the hurdle's weights to this feature, dec-B569 keeps
probit at 0 and 1, and the
2026-08-28 schedule is integer-weight probit's, not the hurdle's.

- (a) 0 and 1 on both parts now, the membership a probit fit takes, read the same under both readings and at predict as
  dec-B580's doors read them; whole numbers when integer-weight probit is built, on probit and the hurdle together; a
  refusal of other values names the remedy, repeating each row its count's number of times. Recommended: the hurdle's zero
  part stays a fit a plain probit call reproduces, so fitting the two parts separately remains the same posterior, as the
  help promises. About 6 agent-hours (S3, subset with it).
- (b) Build integer-weight probit now, for probit and the hurdle together, and take whole-number weights on the hurdle. An
  engine slice in the draws-changing tier (the latent update draws w truncated normals per row; tests/cpp, an exact gate, an
  opus review with mutants), moving the 2026-08-28 schedule and adding probit to dec-B563's families; roughly 31 to 36
  agent-hours, (a)'s 6 and 25 to 30 more, after which frequency-weights' resolution table gains probit.

Evidence that would change it: users of semicontinuous responses with aggregated rows, whose counts (a) takes only by
expanding the rows.

## Constraints

- No engine, bridge, header, state-format or `*.in` file changes; `git diff --stat` on the branch names none.
- Every hurdle call taken today keeps its draws and its methods' outputs: the before/after gate identical, the equivalence
  hurdle scenario bitwise, hurdle-reduction.R bitwise, hurdle-exact.R at its tolerance, test-hurdle.R unchanged. The one
  output allowed to move is a positional ci.level in fitted's third place, which becomes sample (dec-A15).
- Each part stays the standalone probit or gaussian fit at its derived seed, draw for draw, under every new input.
- Each argument evaluated once (dec-B351); both parts see the same rows.
- Out of scope: whole-number or real weights (the weights ruling, dec-B588); a variance forest on the positive part; keepFits = FALSE and
  callbacks; dbarts(), dbartsSpec, xbart and the sampler doors (dec-B14); a sparse design; a two-part formula giving each part
  its own predictors; an engine hurdle (dec-B82).
- Base R calls within DESCRIPTION's R floor; no new dependency.

## Slices

Each slice is a commit series on hurdle-completion; the feature merges once after S6.

### S1: ingestion, the formula interface and the before/after gate (~200 R, ~30 help, ~250 tinytest, ~150 gate)

The hurdle arm ingests once into a data object and derives the parts; the formula and a data object taken, a matrix or
frame x as today; the matrix route's order kept; the zero part's terms without the offset() term (exercised in S4); a data
object's neighbors ignored with dbartsData's warning. Pins: a formula fit equals the matrix fit on the same coded design
draw for draw; factor columns, a level seen only on zero rows (the positive part's test rows code it); row names;
na.action records; the warning beside a data object; the matrix fit of test-hurdle.R unchanged.

The before/after gate, benchmarks/R/hurdle-before-after.R, is written first, before the arm is rewritten. One script, run
as `record <out.rds>` under R_LIBS set to an absolute library path, then `compare <a.rds> <b.rds>`; same machine, two
threads at most. One library is built from this plan's base (main 170635476, shipped mode, `buildInfo()$mode` checked), the
other from the branch. It records, for every call a hurdle fit takes today: a frame x with factor columns; a level with no
training rows (one seen only on zero rows, and one with no rows at all); na.omit, na.exclude and the default with missing
values on zero rows and on positive rows; row names on x; several chains; an unseeded fit under set.seed; n.threads = 2;
the retired spellings the arm restores (sigdf, sigquant, resid.prior, proposal.probs); the residual prior written as
family = hurdle.lognormal(sigma = ), chisq and fixed; keepTrees TRUE and FALSE. For each fit: the fit less its samplers
($zero$fit, $positive$fit, compared through extract(type = "trees") on each part), every draw channel of both parts, and every method's
output: extract and fitted at every type, with ci.level by name and combineChains = FALSE; residuals; extract(type =
"loglik"); predict at newdata (a frame holding the empty levels) at ev, prob, bart and ppd under set.seed, and with no
newdata; pdbart and pd2bart; summary; capture.output of print and summary; family(); plot run on a null device. compare
reports each element by identical() and fails on any difference. Run after each slice and at the end, and by each reviewer
independently from a fresh record on its own libraries.

### S2: the test channel (~150 R, ~45 help, ~380 tinytest)

The stacked positive test, extract and fitted at sample = "test", fitted's formals, keepTrainingFits = FALSE, plot reading
the first n test columns, the test rows' names and weights; man/dbarts.Rd's hurdle sentence. Pins: extract(sample = "test")
at types "prob", "bart" and "ev" equals predict(newdata = test) on a fit kept with its trees, to rounding (the stored draws
against the replayed trees), these types drawing nothing; the training channel unchanged by a test set; the "ppd" test
draws' NA at a 0 test weight; refusal words without a test set. keepTrainingFits = FALSE: its test draws identical() to the
default fit's at one seed, every type ("ppd" under one set.seed); fitted, residuals, extract(type = "loglik"), predict with
no newdata and plot refused in the gaussian fit's words; summary and print answer. plot on a fit with a test set whose first
m rows hold a positive row runs, and its log-scale panel's band equals the one computed from the first n test columns.

### S3: subset and weights (~110 R, ~25 help, ~300 tinytest)

subset on both routes; 0/1 weights; the refusal of other values naming the remedy; a missing weight under na.action
(dec-A89) on both routes; predict's weights formal, at newdata and with no newdata (dec-B470, dec-B482); the 0 row's NA
draws and NaN log-likelihood; partial dependence leaving 0 rows out. Pins: each part equals the standalone fit with the same
subset and weights; a subset's na.action record against lm's on the same rows; residuals at a 0 row; a missing weight on
the matrix route dropped under na.omit and na.exclude, the record against the same fit on the formula route, refused under
the default. With no newdata: under one set.seed, predict(type = "ppd", weights = 1) on a fit with 0 rows gives finite draws
at those rows and, at the weight-1 rows, the draws extract(type = "ppd") gives there; a single 1 and a vector of n ones
agree; a vector of another length is refused in dec-B470's words.

### S4: offsets (~250 R, ~45 help, ~420 tinytest)

The plain and named forms, terms, offset.test, a missing offset under na.action (dec-A89), predict's reserved slot at
newdata and with no newdata (the positive part's test channel shifted), partial dependence. Pins: "prob" moves with a
zero-part offset only and "bart" with a positive-part one only; predict at the training rows as newdata equals fitted where
offsets are expressions; predict with no newdata and the fit's own offset equals fitted; predict with no newdata and an
offset DIFFERENT from the fit's, list(positive = o2) and list(zero = o3), equals predict at the training rows as newdata with
that offset, to rounding, at types "bart", "prob" and "ev", the other part's draws unmoved; a missing value there under
na.pass and na.omit (dec-B482); on a formula fit with an offset() term, predict(fit$zero, newdata) unmoved when the term's
column changes and equal to the hurdle's "prob" at that newdata; refusals naming the part.

### S5: gates (~200 in benchmarks/R)

hurdle-reduction.R gains arms for a formula fit, a test set, a subset, 0/1 weights and each part's offset, each part
bitwise against its standalone fit. The new arms build the standalone inputs by hand, without the package's helpers: the
rows a subset and na.action keep, the 1{y > 0} and log y split, the positive rows' weights and each part's offset are
derived in the script from the data, not by splitHurdleResponse or the arm's own code (a hand-built positive-part formula
fit keeps the same factor levels, c11-formula-levels.R). Today's arm, which takes the split from splitHurdleResponse, stays
as it is. hurdle-exact.R gains an arm with a 0 row in each cell and one with an offset in each part, the expected posterior
derived from the response, weights and offsets alone (the split re-derived, as its header requires). The equivalence hurdle
scenario is not re-recorded.

### S6: records (~50)

Ledger entries for the calls made and the two rulings; this plan's Status and Landing note; TODO; the plans index; NEWS's
hurdle bullet (hurdle is new in 1.0-0, so one bullet against 0.9-34, which had no hurdle family, run in scratch/libs/cran).
The landing tells the maintainer call made 7 (dec-A89 on the matrix route). A TODO pointer is added to the frequency weights
entry: once this lands, frequency-weights' resolution table row for hurdle.lognormal and its refusal 4 (both still read
"weights refused") are edited to the hurdle's 0/1 weights, as the interaction section below says.

## Tests and gates

Per slice, against its own library (`R CMD INSTALL -l <lib> .`, `buildInfo()$mode` checked): the touched test files singly;
the four snapshot files on the reference build; `lintr::lint_package()`, `air format --check .`, check-rc-codoc,
check-win-drift, check-doc-freshness, each on its own exit status; the before/after gate (S1) against the base library.
At the feature's end: the before/after gate again, and by the reviewer from its own records; the full tinytest suite in one
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
11. S4: partial dependence passing the zero part's offset alone - a pdbart "prob" pin with a positive-part offset written
    as an offset() term (no change) and a zero-part one (a change).
12. Any slice: the matrix route's seed derivation or test stacking changed - the equivalence hurdle scenario and the
    before/after gate.
13. S4: with no newdata, the positive part's replacement offset applied through `withTrainingOffset` (its training channel)
    rather than the first n columns of its test channel - the different-offset pin against predict at the training rows.
14. S2: under keepTrainingFits = FALSE, the test columns read from n + 1 as in the default - the bitwise pin against the
    default fit's test draws.
15. S1: the zero part keeping the formula's offset() term - the predict(fit$zero, newdata) pin and mutant 11's pdbart pin.
16. S3: with no newdata, given weights ignored and the fit's own kept - the weights = 1 pin (the 0 rows stay NA).
17. S2: plot selecting over all n + m test columns - the plot pin on a fit whose first m rows hold a positive row.

## RNG class

Neutral. No call taken today moves a draw: today's calls are a matrix or frame x with none of the new inputs, and the
positive part's test with no test set stays the n rows. The equivalence scenario covers one numeric matrix with no missing
values and one chain, hurdle-reduction.R holds by construction and test-hurdle.R's missing-value pins are relational, so the
before/after gate (S1) is what shows the rest: factors, empty levels, every na.action, row names, chains, threads, the
retired spellings and every method's output. fitted's third positional argument changes meaning on hurdle fits only, a
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
man/bart.Rd edits. Either way the feature that lands second states the hurdle's weights.type row; if this lands first, S6
leaves the TODO pointer that frequency-weights' table row for hurdle.lognormal and its refusal 4 need the edit.

dec-B580 names predict(weights = 1) as the door to a 0 row's draw; the ordinal, negative binomial and multinomial predicts
have no weights formal on main or in the frequency-weights plan, which is that plan's to settle.

## Not verified

- That a formula route deriving both parts from one data object reproduces the matrix route's draws for a numeric design;
  S1's first pin shows it or names what differs.
- The exact cost of the n + m wide positive test channel at large n; it is n wide today, and m more columns are what any
  bart fit's test set costs.
- The cost of integer-weight probit in the weights ruling's alternative (b) is scaled from frequency-weights' Polya-Gamma slice and the
  weighted-binary memo, not planned.
- hurdle-exact.R's configuration was read for its header only; S5 confirms a 0 row and an offset fit within its two cells.
- The memory of the stacked channel at large n is the critique's arithmetic, not measured: per draw 2n + |P| + 2m values,
  about 8.6 GB at n = 1e5, |P| = 5e4, m = 1e4 and 4,000 draws, against about 3.5 GB for a gaussian fit with that test set,
  and about 0.6 GB under keepTrainingFits = FALSE.
- brms was read from its source and glmmTMB from its help, neither run.

## The critique and what was amended

A blind critique (scratch/hrd-crit/report.md, probes c01 to c11 on a shipped build of 170635476) found the design sound and
asked for twelve amendments. Each was checked against the code, by a probe where it could be run (a01 to a03 under
scratch/hrd-plan/probes-amend), and all twelve are taken:

1. The positive part's no-newdata offset shifts the first n columns of its test channel; confirmed (a01: `withTrainingOffset`
   refuses an n-long offset on that part; the test-channel shift equals predict at the n rows to 4e-16). Offsets (S4), S4's
   pins, mutant 13.
2. The before/after identical() gate: S1, with the reviewer running it.
3. plot.bartHurdle in the context, door table, S2's pin and mutant 17; man/dbarts.Rd in the help list.
4. keepTrainingFits = FALSE's bitwise test (a02 confirms the premise) and the words its refusals use (c07, a01), mutant 14.
5. dec-B470 cited; the no-newdata weights pin, mutant 16.
6. The data object door (a01: arguments beside a data object are ignored with a warning), call made 6.
7. Open call 1's background: the multinomial matrix shape, the k wording, the pscl offset() difference; glmmTMB taken from
   its help, read for the purpose.
8. dec-A89 on the matrix route (c04b, a03 for 0.9-34), call made 7, told to the maintainer at landing; TODO
   fit-na-action-weights-offset filed with this amendment.
9. The zero part's terms drop the offset() term (a01: a probit fit re-evaluates its term at predict); the predict(fit$zero)
   pin; mutant 11 written with a term, mutant 15.
10. Wording: frame x taken today; predict's 0/1 rule on dec-B580's doors; the window's basis; open call 2's costs in hours;
    the refusal's remedy; dec-B185.
11. hurdle-reduction.R's new arms build their inputs by hand (S5).
12. S6's TODO pointer to frequency-weights' hurdle row and refusal 4.

The critic's five mutants are 13 to 17. The optional saving it raised, building the positive part always with
keepTrainingFits = FALSE, is not taken (call made 8). One point is narrowed rather than taken as written: the critique said
plot's length-n logical could pick wrong columns silently for some m. Recycling over n + m columns keeps the first n
selections right and only adds columns from the first m rows' pattern, so the panel either reads the right columns (no
positive row among the first m) or fails on the length check; the defect is loud, and the fix and its pin stand. On open
call 2 the critique would bring (a) as a confirmation; the call stays open, no ruling settling it.
