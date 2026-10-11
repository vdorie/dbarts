# stan4bart's response range: one mapping, fixed at creation from an initial fit that holds the forest's columns

Status: PROPOSED, 2026-10-10 (dec-B437, dec-B441 to dec-B445). Built by
[stan4bart-creation-mapping.md](../plans/stan4bart-creation-mapping.md). For stan4bart it replaces the
recommendation of [per-chain-responses.md](per-chain-responses.md), whose initial fit left the forest's columns
out. Measured on stan4bart's bartcore branch and dbarts's, arm64 macOS, on a shared machine.

## 1. What the mapping is

stan4bart runs dbarts's BART sampler inside a larger sampler: a parametric mixed model (fixed effects, random
effects, residual sd) drawn by gradient sampling, and a forest. Each sweep the forest is fitted to the response
less the parametric part. dbarts maps what it is given onto [-0.5, 0.5] by two numbers, a lower and an upper end
(the range), and states the leaf prior on that scale.

So the range is two prior parameters. With k = 2, which every continuous stan4bart fit uses unless the user names
another, the prior on the forest's function at any point has its mean at the middle of the range and a standard
deviation of a quarter of its width (width / (2 k)). The parametric part has no intercept, the forest carrying the
level, so the middle of the range is the prior mean of the overall level. A binary response has no range: the
probit scale is fixed.

The width is a regularization dial, and no setting is the true one. At 200 rows a range 10 percent wider cost
1.1 to 1.4 percent of expected-value error in each of five designs and, with many small groups, bought about 4
percent in how well the fit divides a group-level signal between the forest and the random effects. Too narrow,
and the random effects keep signal that belongs to the forest; too wide, and a small sample is overfitted and the
chains mix worse.

Terms used below. "Target width": the width of the response less the generating model's own parametric part.
"Expected-value error": root mean squared error of the fitted mean against the generating mean. "Split error":
root mean squared error of the fitted random effects against the generating ones, both centered. "Today's
method": the branch's, section 2. Differences are paired over seeds on the same data, as percent with a 95
percent interval; under the rule fixed before the fits, "worse" needs the interval to exclude zero on the bad
side and a mean beyond 3 percent.

## 2. Today's method

Each chain is its own sampler. Creation takes the range from the response less the fitted values of an initial
lme4 fit of the model's own terms (lm with the grouping terms as fixed effects where lme4 is not installed).
During warm-up each chain asks dbarts to re-derive the range from the response less its own current parametric
part, at every sweep of the first eighth, every second sweep of the next, and so on, and keeps the last. Of 1000
warm-up sweeps the last re-derivation is at sweep 896: the kept pair equaled that sweep's to 3e-14 in 90 of 90
fits. At the first sweep the installed range was 1.6 to 27 times creation's in every design.

## 3. What was measured

Designs: Friedman's function in ten uniform covariates with noise sd 1, a fixed covariate, and random effects,
at 200, 1000 and 5000 rows, 4 chains of 1000 warm-up and 1000 kept sweeps unless said. "Many small groups" has
groups of 3 to 5, random-intercept sd 4, and two of the forest's covariates constant within group; "control" has
no group-level covariate.

### 3.1 Freezing per chain gives the chains different posteriors

- The chains of one fit kept ranges 1.16 to 1.58 times apart (many small groups, 200 rows, 5 fits); 1.05 to 1.5
  in ordinary designs; 2.8 to 6.6 with 8 groups; with a random slope on a skewed covariate the widest chain kept
  4.5 to 9 creation ranges and chains sat up to 3.4 ranges apart.
- A chain's range decides its answer: across the four chains of a fit the range and the share of the group-level
  signal left in the intercepts correlate at -0.65 to -0.99 in 5 of 5 fits, shares running 0.29 to 0.67 in one.
- More warm-up does not cure it. At 10,000 warm-up sweeps the share's R-hat was 1.25 to 1.44 in 3 of 5 data sets
  under today's method and 1.01 to 1.02 in 5 of 5 under one fixed range.

### 3.2 Fixed at creation from today's initial fit: too narrow

The initial fit has no competitor for a group-level signal, so its random effects take it and the range comes
out narrow: 0.56 to 0.68 of where warm-up settles with many small groups, 0.77 to 0.82 with a confounded
treatment, 0.49 to 0.69 with a random slope, 0.93 to 1.03 in the controls. Against the target width, 0.58 to 0.61 with
many small groups.

- Many small groups (5 seeds; 200 / 1000 / 5000 rows): split error 2.63 / 1.58 / 1.11 against today's 2.14 /
  1.48 / 1.01, paired +22% [-0.2, 44], +7.2% [3.7, 10.8], +9.5% [5.8, 13.2]; coverage of the random effects 24,
  11 and 8 points lower. Unchanged at five times the run length. The expected value is better at 200 rows, -6.1%
  [-9.8, -2.4], and equal above.
- Mechanism: the leaf prior's sd in response units is 0.32 / 0.35 / 0.40 against today's 0.48 / 0.64 / 0.66, and
  the share of the group-level signal left in the intercepts 0.80 / 0.31 / 0.14 against 0.61 / 0.21 / 0.09.
- It is a region: clear at 0.6 of the settled range, marginal at 0.8, gone where the group-level signal is
  halved. Over a later, wider set of 27 design cells it was worse on split error in 14 and on expected-value
  error in none.

### 3.3 Wider fixed ranges overshoot

The response less only the fixed part of the initial fit, and the response's own range, are 1.6 to 2.2 times the
settled range with many small groups and 1.2 to 1.9 times in the controls (702 fits). They repair the split
(1.7 against today's 2.1 at 200 rows) and cost the rest: expected-value error 0.79 and 0.78 against 0.72, its
90 percent interval covering 0.80 against 0.85, effective sample size 165 and 140 against 346.

### 3.4 A range learned in warm-up, per chain or pooled

Learning the range over the second half of warm-up matches today's method at 1000 rows and above and sits halfway
to the narrow fixed range at 200 (split error +10% [-7, 28]). Pooled across chains by the median it is one
mapping, but each chain's draws then depend on the other chains' seeds and on their number; the pooled value
moved a chain's range by up to 48 percent with no warm-up left; a mean was pulled by one stuck chain (three good
chains left with 14 to 16 percent of values outside the range); and chains run as separate processes would have
to stop and exchange values.

### 3.5 A drawn k

With k given dbarts's chi(1.5, 2) hyperprior, two mappings 1.3 to 4.1 times apart in width end with prior sds
within 10 percent of each other in 14 of 19 cells (654 fits in that study): the dependence narrows and does not
go. They differ by 15 percent at 200 rows in both small-group settings and by 29 to 73 percent with a random
slope, where some chains drive k toward zero.

- k mixes badly: 5 to 23 effective draws of 4000, R-hat 1.11 to 2.74; at 5000 + 5000 sweeps, 20 to 73 of 20,000.
- It recovers little or none of the narrow range's loss at 1000 rows and above (split error 1.54 and 1.11
  against today's 1.47 and 1.00): the data choose a prior sd of 3.4 to 3.6 where the split rewards 5.5.
- On whether a drawn k confounds the random-effect variance: it neither biases nor widens the random-intercept
  sd's posterior (interval width 1.65 to 1.74 at 200 rows, 0.65 to 0.67 at 1000), and costs it 2 to 4 times in
  effective sample size where forest and intercepts contest a signal (326 to 78, 456 to 161, 122 to 46) and
  nothing elsewhere. The summed group effect is determined 2 to 4 times better than either part.

### 3.6 An intercept

An intercept in the parametric part, the forest's range centered at zero: the intercept and the forest's level
trade off freely. Chains ended 39 to 80 response units apart in the intercept with their sum unaffected
(expected-value error within 0.036 of today's in 19 of 19 cells), 43 to 66 percent of values lay outside the
range, and k fell to 0.7 to 1.9 with 4 to 5 effective draws.

### 3.7 The initial fit given the forest's columns

Give the initial mixed model the forest's columns as linear fixed effects, and take the range from the response
less that fit's own part only (the model's fixed and predicted random effects, not the added columns'). The
added columns compete for the group-level signal the random effects would otherwise take. A linear term in the
two group-level covariates carries 67 percent of that signal's variance in the main design.

Width over the target, mean of 20 seeds, initial fits only (13,160 of them over 27 designs and 10 rules):

| design, rows | model's own fit | with the forest's columns | lm, no random terms | response's own |
|---|---|---|---|---|
| many small groups, 200 | 0.59 | 0.99 | 1.60 | 1.89 |
| many small groups, 1000 | 0.58 | 0.94 | 1.62 | 1.86 |
| many small groups, 5000 | 0.60 | 0.95 | 1.73 | 1.96 |
| control, 200 | 1.02 | 1.00 | 1.16 | 1.43 |
| confounded treatment, 1000 | 0.78 | 0.91 | 0.97 | 1.13 |
| random slope, skewed, 1000 | 0.81 | 0.80 | 2.24 | 2.61 |
| groups of 10 to 30, 1000 | 0.63 | 0.98 | 1.21 | 1.49 |
| crossed, 1000 | 0.78 | 0.96 | 1.36 | 1.59 |
| weights and offset, 1000 | 0.63 | 0.97 | 1.64 | 1.86 |
| no linear part, 1000 | 0.59 | 0.64 | 1.63 | 1.90 |
| 150 extra covariates, 200 | 0.57 | 1.42 | 1.60 | 1.87 |
| 8 groups, 10 group-level covariates, 200 | 0.62 | 1.32 | 1.35 | 1.63 |

Seed to seed the enriched ratio has sd 0.10 to 0.17 at 200 rows and 0.02 to 0.12 at 1000.

Fits (1,914 in 34 cells; 20 seeds at 200 rows, 10 at 1000, 5 at 5000), against today's method:

- At 1000 rows and above, 18 of 20 cells show no accuracy quantity worse. The two: a group-level signal with no
  linear part (section 7), and a group-level factor entered without indicators, which rule 3 repairs.
- At 200 rows, 4 of 14 cells are worse on expected-value error and none on split error: 150 extra covariates
  (+6.8% [5.5, 8.1]); a group-level signal alone with 50 groups (+4.4% [3.5, 5.3]) and with 8 (+24.2% [16.1,
  32.8]); 8 groups with 10 group-level covariates (+4.5% [1.7, 7.5]).
- Many small groups, 200 rows (20 seeds): expected-value error +2.7% [1.6, 3.9], split error -7.3% [-10.8, -3.8],
  random-intercept sd 4.44 against 4.55 (generating 4), coverage of the expected value 0.842 against 0.860,
  effective sample size 352 against 423. At 1000 rows +0.5% [-2.0, 3.0] on the split and -0.4% [-1.2, 0.5] on
  the expected value. It is a wider, better-centered range than today's 0.85 of the target, not the same answer.
- Per row the expected value's posterior mean moves by a median 0.08 to 0.28 posterior sd, against 0.06 to 0.29
  between two runs of today's method; only many small groups at 200 rows is above its own rerun (0.078 against
  0.058). The rerun was run in 5 of the 34 cells.
- With 1 or 2 chains, and with 100 or 250 warm-up sweeps, it stands to today's method as at the default.
- Noise level: a rerun of today's method at another seed has 2 of 30 intervals excluding zero.

### 3.8 Linear terms or splines for the added columns (dec-B441)

A group-level signal with no linear part earns the forest nothing from linear terms: width 0.63 of the target and
split error +26.2% [16.9, 36.2] at 1000 rows. Three-degree-of-freedom splines per numeric column repair that
(+0.6% [-4.3, 5.8]) and cost at 200 rows with many small groups (+5.9% [4.4, 7.3] expected-value error against
linear's 2.7), with three times the columns.

The in-between was then measured, 240 fits on paired seeds, the rule fixed first: at 1000 rows with the
group-level signal half curved by variance, keep linear if its split error's difference from today's method has
an interval containing zero or a size no larger than the rerun's largest; otherwise splines.

| group-level signal curved, rows (seeds) | rerun | linear | splines |
|---|---|---|---|
| a quarter, 1000 (10) | -0.2 [-2.4, 2.0] | +2.8 [-0.4, 6.0] | -0.6 [-3.6, 2.5] |
| half, 1000 (20) | -1.2 [-3.4, 1.1] | +3.0 [-0.1, 6.2] | +0.2 [-2.2, 2.8] |
| three quarters, 1000 (10) | -0.3 [-4.2, 3.7] | +8.2 [4.7, 11.7] | -1.6 [-5.9, 2.9] |
| all, 1000 (10) | | +26.2 [16.9, 36.2] | +0.6 [-4.3, 5.8] |
| half, 200 (20) | -1.5 [-5.9, 3.0] | -5.9 [-10.8, -0.7] | -23.0 [-27.7, -18.0] |

The rule keeps linear, by a tenth of a point. Linear is inside the rule, not equal to today's method: the share
of the signal left in the intercepts is 3.4 points [2.4, 4.4] higher. At 200 rows half curved, splines cost 6.4%
[4.8, 8.0] of expected-value error, 3.7 points of coverage and 36 percent of effective sample size, linear 1.5%
[0.1, 2.9]. The help says a strongly curved group-level effect is where a user sets the range.

## 4. The design

One linear mixed model by restricted maximum likelihood, from what stan4bart has already parsed from the call:
the model's own fixed columns and random terms as written, plus the forest's columns as linear fixed effects. The
range is the smallest and largest of the response, less any offset, less the model's own part of that fit: its
fixed effects on centered columns plus its predicted random effects. The added columns' contribution is left
out; it is the forest's. Two numbers, computed once, recorded on the fit, the same for every chain.

Building from the parsed pieces and not from a rewritten formula makes `.` inside `bart()`, transformed terms,
`subset`, weights and the offset come out right: the added columns are the columns the forest splits on, on the
rows the model fits.

1. Numeric, logical and ordered columns enter by the value the forest splits on (a logical as 0/1, an ordered
   factor as its level score, a transformed term as transformed). With a logical, an integer-coded category and
   a log term among 15 columns the width was 0.95 to 0.98 of the target; the model's own fit gave 0.70 to 0.73.
2. Aliased columns are dropped, the model's own kept first: a constant, a sum of others, a covariate also in the
   model's fixed part. 3 of 15 were dropped in that design, with no failure in 40 fits.
3. An unordered factor enters as indicators unless every one of its levels sits inside one level of a grouping
   factor; such a factor is left out, since its indicators would take over that random effect. Group-level
   6-level factor, 1000 rows: without indicators width 0.77 and split error +5.7% [1.6, 10.0]; with them 1.00
   and +0.4% [-2.2, 3.1]. The grouping factor itself also inside `bart()`: with its indicators width 1.57 and
   the random effects' effective sample size -22% [-29, -15] for no accuracy gain; left out, 0.93 and no
   difference. A factor finer than the grouping factor was screened by width only (0.95 left out, 1.14 not).
4. The same fit gives the starting values: the forest's first offset is the own part and the first residual sd
   is the fit's. Against starting from the model's own fit under the same range, no accuracy quantity differed
   in 8 cells, warm-ups of 100 and 250 sweeps among them; 6 of the 8 are one design.
5. The range is installed as two numbers before the run and nothing re-derives it. A supplied range and a
   restored fit go through the same step. It needs no change to dbarts.
6. Rows of weight zero are left out of the fit and of the range. Reading the range on the other rows is not
   enough: with a fifth of the rows at weight zero and responses 40 higher there, a mixed-model fit with the
   zero weights left in put the random-intercept sd at 1.9 to 2.1 where the remaining rows give 4.2 to 4.6, and
   the range came out 9 to 31 percent too wide (one design, 5 seeds, responses made deliberately wild).
7. The offset is subtracted once. Today the initial fit's fitted values contain the offset and creation adds it
   again; warm-up hides that, and a range fixed at creation would not.
8. No penalized version. Treating the added coefficients as one more random term puts the width at 0.56 to 0.79
   wherever a few of many columns matter, the region where the narrow range loses.
9. No pilot run. A 100-sweep pilot that sets center and width put the range a full width from the data in 1 of
   20 control fits (95 percent of values outside). One that sets only the width did as well as the enriched fit
   in five hard cells and needs the re-derivation kept, a run before the chains and a seed of its own.
10. Convergence warnings and singular fits do not change the route. Variance parameters moved by a factor of 0.5
    to 2 changed the width by under 2 percent on average and at most 10 percent in one seed (11 cells of 10
    seeds). Singular fits, 23 of 1320, gave the widths of regular ones.
11. No route is chosen by a timer or a size (section 6).

The width as it comes is accepted, with no factor shrinking it toward today's 0.85 at 200 rows: such a factor has
no reason behind it and would undershoot at 1000 rows, where the two already agree.

## 5. Where the fit cannot be computed, and the argument

The fit is undefined with more fixed columns than rows (it stopped in 20 of 20 tries with 400 added columns at
200 rows) and failed numerically in 2 of 1320 fits (1 of 20 at each size with 8 groups and 10 group-level
covariates). Ranges for those cases, on five hard cells against today's method:

| range | width over target | 200 rows, many small groups | 1000 rows and above |
|---|---|---|---|
| enriched fit | 0.91 to 1.00 | expected value +2.7%, split -7.3% | no difference |
| halfway (geometric mean) between the model's own fit and the response's | 0.94 to 1.21 | +4.8%, -8.8%; control +2.8% | no difference |
| lm with the forest's columns, no random terms | 0.96 to 1.68 | +11.0%, -20.9%; effective sample size -62% | accuracy equal; mixing -18 to -29% |
| the response's own | 1.12 to 1.90 | +11.4%, -20.1%; control +5.3% | accuracy equal; mixing -18 to -37% |
| the model's own fit | 0.58 to 1.02 | -5.5%, +20.1% | split +9.5% and +9.9% |

The response's own range, less any offset, is the fallback (dec-B444): it is the definition with the parametric
part taken as zero, it is what dbarts uses alone, and it can always be computed. Halfway is a constant with no
argument and 65 fits as its test. The model's own fit is the one range that loses at every size.

The user's argument, `bart_range` (dec-B445): a pair of numbers, or "mixed" (the default, above), "lm" (the same
fit without the random terms) or "response". They serve a user who knows the range; one whose mixed-model fit
outlasts the sampler (section 6), for whom "lm" costs nothing measured in accuracy at 1000 and 5000 rows and 19
to 24 percent of the random effects' mixing, and 10 percent of expected-value error if used at 200 rows; and one
who expects a strongly curved group-level effect. Splines are not offered, and no size picks a value.

## 6. The fit inside stan4bart (dec-B442, dec-B443)

lme4 is a suggested package, and without it today's initial fit silently becomes another model. With the fit
setting the prior, that would be two models on two machines; importing lme4 would make five more packages hard
requirements. The fit is instead written over Matrix, which stan4bart imports: the profiled restricted
criterion of Bates, Maechler, Bolker and Walker (2015, section 3), from the parsed random-effect structure, by a
sparse Cholesky factor updated at each value of the variance parameters and a bound-constrained optimizer.

A prototype of about 30 lines against lme4 on the same pieces:

- 110 fits over 11 design cells: the own part within 4e-6 of a range's width, in the same time.
- Fifteen further models, among them a correlated intercept and slope, three correlated effects on 15 and on 6
  levels, a slope proportional to its intercept, no group effect at all, nested terms with weights and an
  offset, and crossed factors up to 5000 by 2000 levels: the criterion equals lme4's at three settings of the
  variance parameters to 5e-10 or better (to a constant, the sum of the log weights, under weights); the range's
  ends agree within 2e-5 of the width and the residual sd within 5e-5, the residue being where each optimizer
  stopped (criterion at the prototype's optimum within 2e-4 of lme4's). Five of the fifteen are singular.
- Evaluations: 6 to 400 against lme4's 9 to 136, most where six variance parameters sit on 6 levels.
- On a second draw of seven of those models, one singular fit (the slope proportional to its intercept) stopped
  on the boundary, an intercept variance of zero, 0.017 above lme4's criterion, and did so under each of three
  base-R optimizer settings; the range's ends moved by 2.7e-4 of the width. The other six ended within 5e-5.

Time, one core, seconds; the first three columns are one timing, the last two another, run together:

| rows, grouping | model's own fit | with the forest's columns | one chain, 2000 sweeps | lme4, rerun | inside stan4bart |
|---|---|---|---|---|---|
| 1e5, 25,000 nested groups | 0.2 | 0.4 | 120 | 0.36 | 0.28 |
| 1e5, 5000 by 2000 crossed | 35 | 52 | 65 | 46 | 39 |
| 3e4, 3000 by 3000 crossed | 105 | 111 | 28 | 96 | 76 |
| 1e5, three crossed factors with slopes | over 600 | over 600 | 290 | 9.5 per evaluation | 6.0 per evaluation |

The cost is the factorization of the random effects, not the added columns, wherever the fit is expensive; where
it is cheap many added columns multiply it (10 seconds against 0.3 with 199 columns at 1e5 rows), about 2
percent of a default run. With large crossed factors the initial fit already outlasts the sampler today, and
with three crossed factors and slopes neither routine ends in ten minutes (seven variance parameters took 243
evaluations at a tenth of the size). Fitting the variance parameters on a fifth of the rows saves little (19
seconds against 52): a fifth of the rows still holds most of the levels.

lme4 keeps no run-time use (dec-B443). Two facts about today's package bear on it. stan4bart's copies of lme4's
and reformulas's formula helpers are replaced by those packages' own when the package is installed, not when it
is loaded: a build made beside lme4 and then loaded where lme4 is absent stops at the first call ("could not find
function"), the case of a binary built on one machine and installed on another. And the two builds parse six
formulas (slopes, `||`, nesting, interactions of grouping factors) to identical random-effect structures.

## 7. What does not work

- A group-level signal with no linear part (section 3.8).
- Added columns approaching the number of rows. With 150 extra covariates at 200 rows the width is 1.42 to 1.44
  of the target, expected-value error +6.8% [5.5, 8.1] and its interval's coverage 0.74 against 0.81; with 50,
  1.15 to 1.16.
- Group-level covariates approaching the number of groups. With 8 groups and 10 such covariates the fixed
  columns reproduce the grouping factor's indicators exactly (20 of 20 data sets): the predicted random effects
  are zero whatever the variance parameter, the criterion is flat in it, and the whole between-group variation
  goes to the forest's side of the range (width 1.32 to 1.35). In the sampler's fit the random-intercept sd then
  comes out 1.62 for a generating 4, its interval covering in 7 of 19 fits, where the model's own fit's narrow
  range gives 4.06 and covers in 20 of 20. With 6, 4 and 2 such covariates the widths are 1.22, 1.17 and 1.07.
- Two repairs fail. Leaving out the group-level columns "when they are too many for the groups" turned the
  overshoot into the model's own fit's undershoot (0.61 to 0.64) and its answer moved with the one constant in
  it. Falling back when the enriched fit is singular in a term whose plain fit is not would not fire: in 20
  data sets per cell the enriched fit was singular in 0 with 10 group-level covariates, 1 with 6, 0 with 4 or 2,
  and 0 with 150 or 50 extra covariates; and the response's own range it would fall to is wider still (1.63 and
  1.78) where the overshoot is the fault.
- Eight groups generally. Nothing converged under any method, today's included (R-hat of the random-intercept sd
  1.31 to 1.38), so nothing there can be ranked.
- What the range does not touch: a cluster-level treatment inside the forest is estimated at 3.97 for a
  generating 3 under today's method and under this one alike, covering in 2 of 10.

## 8. Binary responses

No range is derived, and none of the above applies to one. Today a binary fit's starting offset is the initial
fit's fitted probabilities (0.01 to 1.0) where the linear predictor is meant (-2.3 to 3.1 in one check); the
route without lme4 has the same fault. Starting values only: over 30 fits (two designs at 500 rows, 5 seeds, 4
chains of 1000 + 1000) the fit started from glm with the grouping factor as a fixed effect moved the fitted
probability's posterior mean by a median 0.026 to 0.046 posterior sd per row from the fit started from a probit
mixed model, where two runs of the latter differ by 0.025 to 0.040, and its error by +0.1% [-0.6, 0.7] and
+0.1% [-1.5, 1.8]. The glm with a grouping factor as fixed effects is dense in the levels: 353 seconds at 2000
groups of 4 where the mixed model took 0.7, and with groups of 4 its linear predictor is beyond 4 in size on 41
to 44 percent of rows, groups of one outcome being fitted exactly.

## 9. Limits of the evidence

- Nothing is built. The enriched fit, the range's install and the routine are prototypes; the released stan4bart
  was never run, so "today's method" is the branch's.
- About 1,100 quantity-by-cell comparisons were made without correction; single flags near 3 percent are
  unconfirmed. The earlier studies used 5 seeds, where only differences above about 5 percent show.
- The spline question rests on one design family: many small groups, one covariate carrying both parts, one
  curve shape. The factor rules, the zero-weight rule and weights with an offset rest on one design each; the
  factor finer than the grouping factor was screened by width, not fitted.
- Nothing converged in the hard cells of the earlier studies (effective sample size of the expected value 8 to
  70 of 4000 in three of four); large effects survive, claims of equality there are weak.
- Not run: prediction for new groups, non-Gaussian errors, accuracy above 5000 rows, time above 1e5 rows,
  starting values at short warm-up where the overall level wanders.
- The harness removed the whole mean of the own part where the model's centering removes only the fixed
  columns' means; the difference was at most 1.8 percent of a width.
