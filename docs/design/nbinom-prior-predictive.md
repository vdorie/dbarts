# The negative binomial default prior: a prior predictive study

Status: FINDING, 2026-10-10 (dec-B412, dec-B431, dec-B432). A prior predictive study of the default prior
of `family = "nbinom"`, with nine alternative priors on k and fits to five real datasets. The default is
unchanged (dec-B431). The maintainer judged the study flawed and the default improvable; a redesigned
study runs before 1.0-0, after the merge to main (dec-B432).

## Summary

The question (dec-B412): is the default prior of a negative binomial fit, k drawn under
`normal(k = chi(1.5, 2))`, badly off, implying counts and spreads of the log mean no analyst would accept?
Datasets were drawn from the prior over 120 combinations of n, trees, design and mean count, and compared
with other packages' defaults and twelve real datasets; nine alternative priors on k were run through the
same grid and fitted, with the default, to five real datasets. The typical prior dataset is tame and sits
among the real ones; the top tenth does not (one in a hundred has a spread of the log mean of 15, against
a real maximum of 1.4, and 5 to 10 percent at a mean count of 10 hold a count above 10^6, which the
package refuses). By the pre-set rule the default is not badly off, and no alternative fits better beyond
fold noise. The maintainer kept the default (dec-B431), judged the study flawed and the default
improvable, and ruled that a redesigned study runs before 1.0-0, after the merge (dec-B432).

## The prior as built

A row's log mean is a center plus the sum of the trees' leaves. The center is the log of the mean count and
is set by the fit; a user moves it only through an offset. Given k, the leaf scale is set so that the sum
has standard deviation 3 / k at every x, whatever the tree count. k is drawn: k = 2 x chi(1.5), with
median 1.91 and 1st percentile 0.12.

- **Spread of a row's log mean.** Given k, the standard deviation is 3 / k. Over k's prior its median is
  1.57 (a factor of 4.8 on the mean count), its 90th percentile 5.1 (a factor of about 170) and its 99th
  percentile 24 (a factor of about 3e10).
- **Marginal form.** Averaging over k, a row's log mean is the center plus 1.22 times a Student-t variable
  on 1.5 degrees of freedom. The tails are so heavy that the prior mean count is infinite.
- **Dispersion.** The size parameter has a prior on a grid from 1 to 50, with median 8 and never below 1.
  It cannot produce the heavy tail that k does.

The quantity compared below is S, the standard deviation across one dataset's rows of the log mean; the
figures above are for one row given k.

## What was simulated and what it was compared against

The grid has 120 cells: n of 50, 150, 500, 2000 and 10,000; 75 and 200 trees; three designs (3 uniform
columns, 20 uniform columns, 10 columns of which 5 are normal and 5 binary); and mean counts of 0.1, 1, 10
and 1000. Each cell holds 25,000 datasets drawn from the default (5,000 for the alternatives), drawn
through the package.

The references:

- rstanarm's `stan_glm.nb` defaults: normal coefficients scaled by the predictors, exponential(1) on the
  reciprocal dispersion.
- brms's default for the standard deviation of a varying intercept, half-Student-t(3, 0, 2.5), with one
  level per row. Its default on coefficients is flat and has no prior predictive, so it enters only as a
  sensitivity.
- The leaf prior of log-linear BART (Murray 2021), set to its count application's value, a log mean with
  standard deviation 2 at every x (equivalent to k fixed at 1.5), with a beta-prime(5, 3) dispersion.
- Twelve real count datasets (quine, epil, warpbreaks, InsectSprays, discoveries, Traffic, quakes,
  Seatbelts, ships, Insurance, grouseticks, Ornstein).

The pass rule was written before the grid ran. Test 1: in a cell, the default's 99th percentile of S is
more than twice the larger of rstanarm's and Murray's. Test 2: in a cell, more than 1 percent of datasets
have a ratio of largest count to median mean beyond the largest of the references' 99.9th percentiles and
the real datasets' values. The default is badly off if either test fails in more than half the cells (60
of 120). Candidates faced the same tests and a third: a candidate must leave room for every real dataset
below its 95th percentile of S.

## Results

### The bulk

Quantiles of S, as ranges over the three designs (the same at every n, tree count and center):

| prior | 50% | 90% | 99% | 99.9% |
|---|---|---|---|---|
| default | 0.88 - 0.95 | 2.9 - 3.1 | 14 - 15 | 64 - 69 |
| Murray, standard deviation 2 | 1.1 - 1.2 | 1.4 - 1.6 | 1.6 - 2.0 | 1.8 - 2.4 |
| rstanarm | 3.8 - 11 | 6.3 - 13 | 8.4 - 15 | 10 - 17 |
| brms varying intercept | 1.9 | 5.9 | 14 - 15 | 31 - 34 |

The default's median S of 0.93 is the smallest of the four and lies within the real datasets' range of 0.12
to 1.4. rstanarm's median grows with the number of columns (3.8 at 3, 11 at 20).

### The tail

Above the median the default widens fast: 3.1 at the 90th percentile, 15 at the 99th, about 65 at the
99.9th. The consequences for counts, as the share of datasets (percent) in which the largest count exceeds
10^6, the cap the package refuses, and the share whose counts sum past 10,000 at n = 150:

| mean count | count above 10^6, n 50 to 10,000 | sum above 10,000, n = 150 |
|---|---|---|
| 0.1 | 3.2 - 5.2 | 7 |
| 1 | 4.1 - 6.6 | 11 |
| 10 | 5.3 - 8.7 | 23 |
| 1000 | 11 - 18 | 98 |

The refused share is nearly flat in n and rises with the mean count; at a mean count of 10 it is 5 to 10
percent across designs and tree counts.

### Against other priors

At n = 150 and a mean count of 10, a count above 10^6 appears in 6 to 7 percent of the default's datasets,
under 0.1 percent of Murray's, 25 to 99 percent of rstanarm's and 21 percent of the brms varying-intercept
default's. The default is far tamer than rstanarm's in the bulk and comparable to brms's in the tail; only
Murray's calibration is thin-tailed.

The pass rule, as written, returns "not badly off". Test 1 fails in 0 of 120 cells; the default's 99th
percentile of S is between 0.86 and 1.8 times the larger reference's. Test 2 fails in 2 cells as written
and in 40 after the rule's error about the count generator is corrected (see the limits); all 40 are the
3-column design, where 1.6 to 2.2 percent of datasets are beyond reach against the 1 percent limit. Both
are under the 60 needed.

### Real datasets

S of each real dataset, from a negative binomial fit by maximum likelihood, and its percentile under each
prior at that dataset's design, n and center:

| data | S | percentile, default | percentile, Murray | percentile, rstanarm |
|---|---|---|---|---|
| quine | 0.41 | 7.6 | 0.06 | 0 |
| epil | 0.79 | 41 | 9.4 | 0.04 |
| warpbreaks | 0.23 | 2.5 | 0.22 | 0 |
| InsectSprays | 0.82 | 51 | 28 | 0.02 |
| discoveries | 0.17 | 0.18 | 0 | 5.4 |
| Traffic | 0.12 | 0.02 | 0 | 0 |
| quakes | 0.49 | 21 | 0.02 | 0 |
| Seatbelts | 0.15 | 0.14 | 0 | 0 |
| ships | 0.55 | 21 | 0.68 | 0 |
| Insurance | 0.31 | 2.2 | 0 | 0 |
| grouseticks | 1.4 | 67 | 66 | 4.7 |
| Ornstein | 0.79 | 58 | 37 | 0 |

No real dataset falls above the 67th percentile of the default. Under rstanarm's default 10 of the 12 fall
at or below the 0.04th percentile; the real data are all far tamer than that prior. Several fall in the
lowest few percent of the default as well (Traffic, Seatbelts, discoveries), so the default is also wide
for those.

## The alternatives and the fits

Each alternative changes only the prior on k; the dispersion prior and the leaf scale stay. What each sets:

- k fixed at 3, 2 or 1.5 (a standard deviation of the log mean of 1, 1.5 and 2).
- k drawn as chi(3, 1.25), chi(5, 0.9) or chi(10, 0.62): more degrees of freedom, with the scale chosen to
  hold k's median at 1.9, as in the default.
- k drawn as chi(1.5, 3) or chi(1.5, 4): the default's shape at two-thirds and one-half of its width.
- k drawn as chi(3, 2): the default's scale with 3 degrees of freedom, which thins only the small-k side.

### Prior predictive results

Quantiles of S, and the share (percent) of datasets with a count above 10^6 at a mean count of 10, as a
range over cells:

| prior on k | S 50% | S 90% | S 99% | count above 10^6 |
|---|---|---|---|---|
| chi(1.5, 2), default | 0.93 | 3.1 | 15 | 5.1 - 10 |
| 3 fixed | 0.58 | 0.74 | 0.91 | 0 |
| 2 fixed | 0.87 | 1.1 | 1.4 | 0 |
| 1.5 fixed | 1.2 | 1.5 | 1.8 | 0 - 0.02 |
| chi(3, 1.25) | 0.92 | 1.9 | 4.4 | 1 - 3.4 |
| chi(5, 0.9) | 0.94 | 1.6 | 2.8 | 0.16 - 1.4 |
| chi(10, 0.62) | 0.93 | 1.4 | 2.0 | 0 - 0.14 |
| chi(1.5, 3) | 0.62 | 2.1 | 9.3 | 2.8 - 6.1 |
| chi(1.5, 4) | 0.46 | 1.5 | 7.0 | 1.8 - 3.8 |
| chi(3, 2) | 0.58 | 1.2 | 2.7 | 0.2 - 0.9 |

Eight of the nine pass the rule. k fixed at 3 fails its third test: it leaves too little room for epil,
grouseticks and Ornstein, whose S sits above its 95th percentile. The drawn priors with more degrees of
freedom keep the default's median S and cut the 99th percentile to 2.0 to 4.4. Grouseticks, the real
dataset with the widest spread, falls at the 67th percentile under the default, 77th under chi(5, 0.9),
89th under k fixed at 2 and 90th under chi(3, 2); a thinner prior leaves less room above it.

### Held-out fits

The default and five alternatives were fitted to five real datasets (Traffic, Insurance with an offset,
quine, epil, grouseticks; 1,033 rows), 4 chains of 500 burn-in and 500 kept draws, scored by 5-fold
held-out log score, five repetitions. Entries are the candidate's score minus the default's, mean over
repetitions (positive means the candidate fits better):

| data | chi(5, 0.9) | chi(10, 0.62) | k = 2 fixed |
|---|---|---|---|
| Traffic | +0.6 | +1.1 | +2.2 |
| Insurance | -0.8 | -1.5 | -3.5 |
| quine | +0.2 | +0.2 | +0.3 |
| epil | +0.7 | +0.5 | +1.0 |
| grouseticks | -0.9 | -4.1 | -11.4 |
| sum | -0.2 | -3.8 | -11.4 |

Two further drawn priors differ from the default by -2.2 (chi(3, 2)) and -0.4 (chi(3, 1.25)) in sum,
neither beyond noise.

Fold-to-fold standard deviation of a dataset's score is 5 to 10 units, so most differences are inside
noise. The pattern is monotone: the tighter k is held near 2, the more is gained on Traffic and the more
lost on Insurance and grouseticks. k fixed at 2 and chi(10, 0.62) are reliably worse on Insurance and
grouseticks (the same sign in all five repetitions). No alternative is better than the default beyond
noise.

### Small n

At n = 30 (Traffic, grouseticks, and simulated wide data with a true S of 3), the
default's tail does not reach the posterior: the mean posterior probability that the forest's standard
deviation exceeds 5 is at most 0.11 percent in every setting, against about 10 percent in the prior.
Thin tails cost on wide data: k fixed at 2 gives 90 percent interval coverage of 0.65 against 0.92 for the
default, and an RMSE of the log mean of 1.31 against 0.94, in all five replicates.

## The decision

The maintainer, choosing between keeping `normal(k = chi(1.5, 2))` and changing counts to chi(5, 0.9)
(dec-B431): "Keep chi(1.5, 2), but I would note a few things: the default already changes
between different families (i.e. k = 2 for Gaussian), it seems as if there is room to improve the default,
and it seems as if the study as flawed." The default is unchanged. The three notes stand beside it: a
default that differs by family is no cost in itself; the default may be improvable; and the study does not
settle that.

On whether a redesigned study runs before 1.0-0 (dec-B432): "Yes, but we can do it post-main merge." It
runs before 1.0-0 and after the merge to main; any change of default it supports comes back to the
maintainer.

## Limits of this study

These are the reasons the maintainer does not take the result as settling the default.

- **The yardstick was other packages' defaults.** The pass rule asked whether the default is much wider
  than rstanarm's and Murray's. Those are themselves loose: rstanarm's median S is 4 to 11 and the
  varying-intercept default's 99th percentile equals the default's. Passing against them says little about
  whether the prior is reasonable in absolute terms. No criterion was fixed in absolute terms beforehand.
- **The pass rule contained an error.** It counted any count above 2^31 - 1 as unfittable, on the belief
  that R's count generator fails there. It does not. Corrected, test 2 fails in 40 of 120 cells instead of
  2; both readings are under the 60 that make the default badly off, and both are reported above.
- **The singled-out alternative was chosen after seeing fits.** The rule chose two candidates for the fits
  (k fixed at 2 and chi(10, 0.62)). The other three, chi(5, 0.9) among them, were fitted after those, and
  chi(3, 2) was devised to fill a gap the first fits showed. The alternative put forward, chi(5, 0.9), was
  picked after seeing all of these fits.
- **Few and small datasets.** Fits used five textbook datasets (1,033 rows in all) and placement twelve.
  None has the wide spreads real count data can have; the widest has S of 1.4.
- **k did not mix.** Under every drawn-k prior at the chain lengths used, k's effective sample size is
  about 6 to 24 of 2,000 and its R-hat exceeds 1.1 in 6 to 10 of 10 seeds on every dataset. Fit comparisons rest
  on noisy posterior means and on held-out scores whose differences are mostly within fold noise.
- **Only the prior on k varied.** The leaf scale (3 / k) and the dispersion prior stayed fixed, so the
  study cannot say whether the improvable part is k's prior.
- **A narrow grid.** No offsets and independent columns only; the default tree prior throughout; the 99.9th
  percentiles rest on 25 of 25,000 datasets a cell; the small-n settings are five replicates each, and the
  wide simulation is one function.

## What a redesigned study should do

- Fix the criterion first and in absolute terms: what a prior dataset should look like, without reference
  to other packages' defaults.
- Use count data with wide spreads, simulated and real, so a thinner prior can be shown to cost something
  or nothing.
- Run chains long enough for k to mix, or fit with k integrated out, so fit comparisons are not noise.
- Vary more than the prior on k: the leaf scale, the dispersion prior and the tree prior.
- Include offsets and correlated columns in the grid.
- Choose any candidate before seeing its fits, and write the rule for choosing it down first.

## Evidence

The study's scripts and outputs are kept outside the tree. They hold the simulation grid, the pass rule as
written before it ran, the candidate list, the placement of the twelve real datasets, and the per-fit
results behind every table above.
