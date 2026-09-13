# Is chi(1.5, 2) still the right binary k prior?

agent: measurement
rng: neutral (nothing in the package changes; this is a measurement)
budget: one benchmark script, this note

## Goal

The binary end-node prior on k moved from chi(1.25, Inf) to chi(1.5, 2)
on the strength of a study that held the degrees of freedom at 1.5,
varied only the scale, and ran four simulated data-generating
processes. The degrees of freedom were never studied, and four
simulated processes are a narrow base for a default every probit fit
rides. This note re-opens the question with both parameters varied and
a much wider case set, and says whether the shipped default should
move. The harness is benchmarks/R/binary-hyperprior.R; it is checked in
and nothing in it changes a package default.

## The grid

Twenty-three priors. Twenty hyperpriors, chi(df, scale) for df in
{1, 1.25, 1.5, 2, 3} crossed with scale in {1, 2, 5, Inf}, and three
fixed values of k, 1, 2 and 3. Fixed k = 2 is what BayesTree used, what
dbarts still uses for continuous responses, and the value chi(1.5, 2)
is centered near.

A hundred and sixty-eight cases.

- A hundred and sixty-two simulated cells: six processes crossed with
  n in {100, 500, 2000}, p in {5, 20, 50} and a base rate in
  {0.05, 0.2, 0.5}, eight repetitions each. The processes are the
  Friedman surface; a sparse linear one; a smooth additive one at a
  high signal-to-noise ratio and the same surface at a low one (the
  "strong" and "weak" pair); a near-separable linear one, scaled so
  that almost every true probability sits within a hundredth of zero
  or one; and one built from two-way interactions with no main effects
  worth speaking of. Past the fifth column every predictor is noise, so
  p = 20 and p = 50 are mostly noise. The base rate is set by an
  intercept solved once on a large fixed pool, so a process is a fixed
  truth rather than something that moves with the training draw. Each
  repetition draws its own training set and its own thousand held-out
  rows, and the true probability of every held-out row is known.
- Six real datasets, scored by sixty repeated eighty-twenty splits:
  Pima (532 rows, 7 predictors, 33 percent positive), biopsy (683, 9,
  35), infert (248, 5, 33), kyphosis (81, 3, 21), pbc dichotomised at
  death (276, 16, 40) and birthwt (189, 8, 31). Two of them carry
  factor predictors. There is no ground truth here, so only the two
  outcome-based scores are available.

Within a case and a repetition every arm sees the same data and starts
from the same seed, so the arms are paired and a difference between
them is a difference in the prior. Each fit uses the package defaults
except the chain count: seventy-five trees, five hundred draws kept
after five hundred discarded, one chain, one thread. That is 38,088
fits at about 0.42 seconds each, and the prior buys none of that time:
0.416 seconds a fit at scale 2 against 0.427 at scale 1.

The scores, all on held-out rows: log score and Brier against the
held-out outcome everywhere; on the simulated cells, additionally the
coverage and mean width of the 90 percent posterior interval for the
true probability, and the root mean squared error of the posterior mean
probability. Standard errors below are over the 1,296 simulated
(cell, repetition) pairs or the 360 real (dataset, split) pairs. The
paired standard error is the one to read: the unpaired one is dominated
by how far apart the cases are from each other, which no prior affects.

## What the scores say

Each row is one score family. "Best" is the arm with the lowest
average; the last column is that arm's paired difference from the
incumbent, which is what says whether the gap is real. Interval
coverage is scored as the distance from the nominal 0.90, so that
over-coverage is not free, with the incumbent's mean width alongside
for what the distance costs.

| score | best arm | best (SE) | chi(1.5, 2) (SE) | paired difference (SE) |
|---|---|---|---|---|
| simulated log score | chi(3, 1) | 0.3020 (0.0047) | 0.3027 (0.0047) | -0.00068 (0.00028) |
| simulated Brier | chi(3, 1) | 0.09521 (0.00172) | 0.09542 (0.00172) | -0.00021 (0.00008) |
| simulated probability RMSE | chi(3, 2) | 0.1199 (0.0018) | 0.1203 (0.0017) | -0.00036 (0.00022) |
| simulated coverage, distance from 0.90 | chi(1.25, 1) | 0.1362 (0.0042) | 0.1533 (0.0047) | -0.0170 (0.0026) |
| simulated interval width, same two arms | chi(1.25, 1) | 0.3433 (0.0058) | 0.3276 (0.0054) | +0.0157 (0.0011) |
| real log score | k = 2 | 0.4290 (0.0095) | 0.4325 (0.0095) | -0.00358 (0.00068) |
| real Brier | k = 2 | 0.1406 (0.0034) | 0.1416 (0.0034) | -0.00106 (0.00021) |

Read down the first three rows first: on everything that scores a point
prediction, the fifteen finite-scale hyperprior arms are one
undifferentiated block. The best of them beats chi(1.5, 2) by two and a
half standard errors on log score and by a quantity - seven
ten-thousandths of a nat - that no user would notice, and the ordering
within the block changes from score to score. What separates cleanly is
everything outside that block. Every improper-scale arm is worse than
the incumbent on log score - by 0.0028 at df = 3, by less than a
standard error at df = 1.25 - and worse on coverage distance at every
df by four standard errors or more. Fixed k costs 0.0036 in log score
at k = 2, 0.0047 at k = 1 and 0.0125 at k = 3, each many standard
errors.

The fixed-k penalty is a small-sample penalty and nothing else. Against
chi(1.5, 2), fixed k = 2 is worse by 0.0102 in log score at n = 100
(SE 0.00093) and by 0.0033 at n = 500, and better by 0.0027 at
n = 2000 (SE 0.00038). Whatever else the hyperprior is buying, at two
thousand rows it is not point-predictive accuracy.

The real datasets say the same thing more sharply, and they are the
most uncomfortable result here. Fixed k = 2 is the best arm on both
outcome scores, better than chi(1.5, 2) on all six datasets and by more
than two standard errors on three of them (infert, kyphosis, pbc). The
mechanism is visible in the sampled k: under chi(1.5, 2) its posterior
median runs from 1.08 on biopsy to 2.85 on birthwt, and the datasets
where the hyperprior loses are the ones where it settles below 2. On
these six datasets the hyperprior shrinks too little, at a cost of
about four thousandths of a nat, under one percent of the score.

Coverage is where the arms genuinely differ, and also where all of them
are in trouble. Against a nominal 0.90, chi(1.5, 2) averages 0.80,
fixed k = 2 averages 0.66 and fixed k = 3 averages 0.54. The shortfall
is not spread evenly over the processes: for chi(1.5, 2) it is 0.97 in
the weak-signal cells and 0.92 in the linear ones, against 0.72 in the
strong-signal cells and 0.54 in the near-separable ones. Near-separable
truth is where a BART probability interval is too narrow, and no prior
in this grid repairs it - the best-covering arm reaches 0.55 there.

One more fact governs how to read the coverage column. Rank the twenty
hyperprior arms by mean coverage and by mean width and the two
orderings agree almost exactly (rank correlation 0.99). Within the
hyperprior family, covering better is being wider, and there is no
interior optimum: the grid's widest prior is its best-covering prior,
and a grid extended to smaller scales would keep going. The fixed arms
are what break that relation and are therefore the argument for a
hyperprior at all: fixed k = 1 has wider intervals than chi(1.25, 1)
(0.368 against 0.343) and worse coverage (0.793 against 0.822),
because it is wide everywhere instead of wide where the data ask.

## Where each prior is worst

Worst-cell regret: average each arm over a cell's repetitions, subtract
the best arm's average in that cell, and report the largest shortfall
any of the 168 cases imposes. This is the robustness column, and it is
the one the July study decided on.

| prior | mean regret | worst regret | worst case |
|---|---|---|---|
| chi(3, 1) | 0.0047 | 0.0165 | weak, n = 100, p = 5, rate 0.5 |
| chi(1.25, 2) | 0.0054 | 0.0193 | weak, n = 100, p = 5, rate 0.05 |
| chi(1, 5) | 0.0055 | 0.0198 | strong, n = 100, p = 5, rate 0.2 |
| chi(1, 2) | 0.0052 | 0.0200 | interaction, n = 2000, p = 50, rate 0.5 |
| chi(1.5, 2) | 0.0054 | 0.0214 | separable, n = 100, p = 50, rate 0.2 |
| chi(2, 2) | 0.0051 | 0.0217 | separable, n = 100, p = 50, rate 0.2 |
| chi(1.5, 1) | 0.0050 | 0.0229 | weak, n = 100, p = 5, rate 0.5 |
| chi(1.25, 5) | 0.0053 | 0.0251 | strong, n = 100, p = 5, rate 0.2 |
| chi(3, 2) | 0.0055 | 0.0274 | separable, n = 100, p = 50, rate 0.2 |
| chi(2, 5) | 0.0057 | 0.0277 | interaction, n = 100, p = 20, rate 0.5 |
| chi(1.5, 5) | 0.0057 | 0.0322 | friedman, n = 100, p = 50, rate 0.5 |
| chi(1, Inf) | 0.0061 | 0.0337 | friedman, n = 100, p = 50, rate 0.5 |
| chi(2, 1) | 0.0053 | 0.0352 | weak, n = 100, p = 5, rate 0.05 |
| chi(1, 1) | 0.0061 | 0.0352 | weak, n = 100, p = 5, rate 0.05 |
| chi(1.25, 1) | 0.0054 | 0.0352 | weak, n = 100, p = 5, rate 0.05 |
| chi(3, 5) | 0.0068 | 0.0453 | friedman, n = 100, p = 20, rate 0.5 |
| chi(1.25, Inf) | 0.0059 | 0.0466 | friedman, n = 100, p = 50, rate 0.5 |
| chi(1.5, Inf) | 0.0064 | 0.0499 | friedman, n = 100, p = 50, rate 0.5 |
| k = 1 | 0.0102 | 0.0665 | weak, n = 100, p = 5, rate 0.5 |
| k = 2 | 0.0087 | 0.0718 | separable, n = 100, p = 50, rate 0.5 |
| chi(3, Inf) | 0.0083 | 0.0728 | friedman, n = 100, p = 50, rate 0.5 |
| chi(2, Inf) | 0.0070 | 0.0941 | friedman, n = 100, p = 50, rate 0.5 |
| k = 3 | 0.0174 | 0.1341 | separable, n = 100, p = 50, rate 0.5 |

The column separates three groups and not more. Every finite-scale
hyperprior's worst case costs under 0.036; every improper-scale arm
except df = 1 costs more than that, up to 0.094; every fixed k costs
0.067 or more, up to 0.134. The ordering inside the finite-scale group
should not be read: it is a maximum over 168 cases at eight
repetitions, the noisiest statistic in this note, and the arms it puts
first and last are indistinguishable everywhere else.

Two other things the column shows. Every worst case is at n = 100 -
which is what makes the small-sample finding above the whole of the
robustness story - and they cluster in the near-separable, weak-signal
and high-noise-predictor corners. And the worst cases of fixed k = 2
and k = 3 are both the near-separable cell, where forcing a large k
over-shrinks a signal that is nearly deterministic. That is the same
failure the July study reported at the same kind of cell.

## What the degrees of freedom do

At a finite scale, the degrees of freedom and the scale are close to
redundant. Both move the same quantity, the prior median of k, which is
the scale times the square root of the median of a chi-squared on df
degrees of freedom: 0.67 at df = 1, 0.82 at 1.25, 0.95 at 1.5, 1.18 at
2 and 1.54 at 3. The fit responds to the product and not to the two
separately. Holding scale at 2 and raising df from 1 to 3 raises the
average sampled k from 1.50 to 1.70, narrows the average interval from
0.332 to 0.320 and lowers coverage from 0.809 to 0.784, while log score
moves by 0.0003 and Brier by 0.00007, neither of them a standard error.
chi(3, 2) and chi(1.5, 5) have nearly the same sampled k (1.70 and
1.78), the same coverage to within half a point, the same width to
within one percent and identical log scores, and they differ in both
parameters. No df in the grid is distinguishable from 1.5 on any
point-prediction score at any scale.

What df does on its own is set how much prior mass sits near zero. A
chi with one degree of freedom has its mode at zero and lets k go small
when the data ask; at df = 3 the mode is well away from zero and k
cannot. That shows up where small k is right: in the near-separable
cells at n = 100, chi(1.5, 1) samples a median k of 0.92 and chi(3, 2)
of 1.27, and the low-df arms cover those cells better. It is a real
effect and a small one, worth a percent or two of coverage at a
comparable cost in width.

At the improper scale df is the entire prior, and there it matters a
great deal. The average sampled k is 5.3 at df = 1, 9.9 at 1.5, 65 at
2 and 305 at 3, with df = 1.25 out of order at 25 because a handful of
runaway chains dominate the average and the ordering among runaways is
noise. Worst-case log-score regret rises from 0.034 at df = 1 to 0.073
at df = 3. The 1e6 cap that the manual calls a backstop that only
engages at the improper scale did engage: three of the 38,088 fits hit
it, all of them chi(3, Inf), all in weak-signal cells. The manual's
description is exactly right and the cap is not decorative.

So the degrees of freedom are a fine adjustment of the same thing the
scale sets, with no value in the grid separable from another on
predictive accuracy, plus a genuine but second-order effect on how far
down k can go; and at an improper scale they are what decides how badly
k runs away, larger being worse. Neither finding argues for moving off
1.5.

## Recommendation

Keep chi(1.5, 2). The evidence that decides it is the shape of the
first table rather than any single number in it: across a case set
three and a half times the size of the one that set the default, and
with the degrees of
freedom varied for the first time, every finite-scale hyperprior scores
the same on everything that scores a point prediction, the incumbent
included, and the only score that separates them - interval coverage -
separates them monotonically in interval width, so it names no interior
optimum and would send the default to the smallest scale on the grid
and then off it. A default cannot be chosen on a criterion with no
stopping point, and the two criteria that do have one pull in opposite
directions and are both small: coverage would move the scale down to 1,
worth 1.7 points of coverage distance for five percent wider intervals,
while the real datasets and every simulated cell at n = 2000 would move
it up, since fixed k = 2 beats the hyperprior on both outcome scores
there by shrinking more than the sampled k does. chi(1.5, 2) sits
between them, and its prior median of 1.9 is the field's k = 2 sampled,
which is the property that made it defensible in the first place. What
the wider evidence does settle, and settles more firmly than the July
study could, is the two negatives: the improper scale is worse than
the finite ones at every df - decisively on coverage and on worst-case
regret, and on log score too once df reaches 1.5 - and no fixed k is
tolerable across the size range, since fixed k = 2 pays 0.010 in log
score at n = 100 and 0.09 to 0.17 in coverage distance at every size.
The hyperprior earns its place; which finite scale and which df it uses
does not matter nearly as much as the default's history suggests.

## What this cannot settle

- Absolute coverage. No arm reaches nominal, and the deficit is
  concentrated in near-separable and strong-signal truth, where the
  best arm covers 0.55. Whether that is the leaf prior, the number of
  trees, the probit link or the forest itself is not separated here,
  and it is a larger effect than anything the prior choice controls.
- Whether the sampled k is systematically too small on real data. Six
  datasets and one dichotomised survival outcome, with no ground truth
  and only outcome-based scores, is enough to notice the pattern and
  not enough to act on. It is the one result that would change the
  recommendation if it held up over a few dozen datasets, and the
  change it would argue for is a larger scale, not a smaller one.
- The ordering among finite-scale arms, on any score. Eight
  repetitions resolve the block-to-block differences reported above
  and do not resolve the within-block ones; more repetitions would
  narrow the paired standard errors but the differences they would
  resolve are already smaller than a user could act on.
- Anything off the fitted configuration: one chain of five hundred
  draws, seventy-five trees, package defaults otherwise. The
  tree count in particular interacts with k by construction, and this
  note holds it fixed.
- The logistic and weighted binary paths. Only the probit family was
  fitted.

## Re-running

    Rscript benchmarks/R/binary-hyperprior.R blocks
    Rscript benchmarks/R/binary-hyperprior.R sim:friedman:500 <outdir>
    Rscript benchmarks/R/binary-hyperprior.R real:biopsy <outdir>
    Rscript benchmarks/R/binary-hyperprior.R summarize <outdir>

One block is one invocation and writes one rds, so the grid splits
across a session; `all` runs every block in sequence. The full run
above is twenty-four blocks and about ninety minutes on four cores. The
tables in this note are the `summarize` output, rounded.
