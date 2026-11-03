# The negative binomial default prior: a study with its criterion fixed first

Status: FINDING, 2026-10-31, corrected 2026-11-03 (dec-B432, dec-B502, dec-B556, dec-B558). The redesigned study of
the default prior of `family = "nbinom"`, run as pre-registered. Its first look, a reading added after it and a fresh
block on new replications all leave the rule keeping the default, k drawn under chi(1.5, 2). The maintainer kept it
for now with the question open (dec-B558); nothing is built.

## Summary

The rule keeps the default, and that default fails the prior check the maintainer ruled. The maintainer kept the
default for now and left the question open.

On a negative binomial fit the forest models the log mean, and k sets how far that log mean can spread across rows: a
small k allows a wide spread. The default draws k from chi(1.5, 2). The rule was fixed before any fit was scored. A
check on the prior alone asks, among other things, that at most 1 percent of the datasets the prior generates have a
spread of the log mean above 4 (Part A). Fits then read intervals on known truths, and real outcomes (Part B). A
candidate that passes the prior check replaces a default that fails it only if it loses nothing on the fits.

The default fails the prior check: 4.0 to 8.4 percent of its datasets pass a spread of 4. Five registered candidates
pass, and all five lose on the fits. A cell below is one dataset or design at one spread.

- **Tighter chi priors, chi(4, 1.03) and chi(5, 0.9),** have about the default's median with less room near zero.
  They fail the check that intervals on narrow data grow by at most 10 percent: on null data, where the forest has
  nothing to find, intervals are 21 and 26 percent wider. They also lose coverage of group means on one and two small
  samples of the `epil` data.
- **A floor on k, at 0.6 and at 0.75,** leaves narrow data exactly as it is. These fail on the interval score of group
  means, on average 1.25 and 1.39 times the default's against a limit of 1.05. That average comes from 4 of 68 cells,
  all widened or small-sample versions of the `epil` and `mite` data, at 2.0 to 5.4 times the default's; the median
  cell's ratio is 1.000 and no other cell is above 1.05. They also lose coverage of group means on two and four cells.
- **A half-normal prior with median 1.3** under-covers one narrow cell, the Traffic data with a BART-fitted truth (0.60
  against the default's 0.78), and loses coverage of group means on two `epil` cells.

The coverage of group means the five lose on `epil` is mostly on groups whose true mean count is near zero. Of the 42
group readings the floor at 0.6 loses there, 35 have a true mean below 0.01, with a median of 0.0003, and on 40 the
truth lies below the interval: a floor on k stops the fit reaching that far down. The plan scores rows below 0.01
apart for the log mean at a row and for a new count, but not for group means. With them set apart for group means
too, a line drawn only after seeing the scores, the floor at 0.6 would no longer fail; [Sensitivity](#sensitivity)
reports it and does not offer it as an option.

A change would buy the prior check and, on long runs, lower interval scores on the log mean and on rate ratios. About
70 percent of the log mean's fall comes from one design, 60 rows at a spread of 3, where the default's 90 percent
intervals cover only 0.54. Without that design the log mean's ratio is 0.82 to 0.92, from the 10 to 14 of 68 cells
below 0.95, and the median cell's is about 1.

A sixth prior, chi(2.75, 2), found afterward by a search along a raised scale, was then tested on 359 new
replications, with both readings stated in advance. It passes on all replications together and in the default-length
check. On the new replications alone two statistics are not shown under their limits: the interval score of group
means, 1.02 against a limit of 1.05 (z 1.37 against 1.96), the same few small-sample `epil` cells losing again; and the
coverage gap of rate ratios, +0.007 against 0.02 (z 1.73), from one replication of one designed cell. So the rule
still keeps the default.

Read beside the rule, these results are weaker than they look, in both directions:

- **That no cell is worse by the bar means little.** At 5 replications a cell, a loss of the same size in each
  replication it touches is flagged only if it touches at least 4 of the 5.
- **What chi(2.75, 2) buys is mostly a long-run figure.** At `bart()`'s default lengths, which users run, its interval
  scores are 0.97 of the default's on the log mean and 0.94 on rate ratios, against 0.85 and 0.85 on long runs, and
  1.03 on group means.
- **The statistic that kept the default rests on a few synthetic cells.** Three small, widened cells carry 98 percent
  of its noise. Without `epil`'s cells it reads 0.95, and without the one design that pulls the other way 1.05; setting
  near-zero groups apart moves it non-monotonically. A block of this size had only a 55 to 65 percent chance of
  passing it even if the two priors were equal.
- **At the spreads of real data the two priors barely differ.** The median cell's interval score ratio is 0.99 to
  1.00 in each class, coverage moves by at most 0.006, and a posterior mean moves less with the prior than between a
  default-length and a long run under the same prior.
- **The repair the prior check asks for moves the bulk of the prior,** not only its tail: the median of k goes from
  1.9 to 2.9.

Read after the fact and deciding nothing, the point estimates and effects users report do not separate the priors
either. For the log mean, group means and rate ratios (per row, within groups and on average), on real data at its
own spread, widened data and designed cells, RMSE ratios run 0.96 to 1.05 and bias, small to begin with, moves by
under 10 percent either way. Only on the small, widened samples does chi(2.75, 2) gain, on the log mean and on
effects, and more on long runs than at default length.

R's truth on seven designed cells was first stored as zero by an error in the study's code, though their counts carry
the effect; every rate-ratio reading here is scored against the corrected truth, and no verdict depends on it.

Fixed k that passes the prior check, 1 to 1.2, fails the narrow-data guard by far, at about 1.9 times the default's
interval score on null data.

The first look used 94.5 of the 110 core-hours allowed, an added reading brought that to 109.9, and a fresh block, on a
stop the maintainer raised, to 137.0 of 137. The rule funds no further runs. A study that reopens the question would
need datasets people meet, chosen by a stated rule and read at their own spreads; estimands natural to each dataset,
with RMSE, bias and effects scored; a reading of what users run, decided first; enough replications to decide; and a
criterion argued as others argue for a default prior ([The outcome, and what a reopened study would need](#the-outcome-and-what-a-reopened-study-would-need)).

## The question and the rule

The question: is there a better default prior on the forest's spread for negative binomial fits,
judged by a criterion fixed before any fit is scored? The criterion has two parts. Part A reads the prior alone.
Part B reads fits: coverage of 90 percent intervals on truths with known values, the interval score, and
held-out log score. A candidate that passes A replaces a default that fails A only if it loses nothing on B. The
rule, its bars and its bounds are as pre-registered; the records are in [Appendix: evidence](#appendix-evidence).

Part B's estimands: **F**, the log mean at a row; **G**, the mean count over a covariate-defined group; **R**,
rate ratios of a binary covariate; **P**, a new count at a held-out row. A class's average is the mean of four
spread bands' averages, cells equal within a band: narrow, middle, wide and small (small-sample wide).

The long reading is each candidate's posterior from runs long enough for k to mix, on a scratch build carrying
the collapsed scale move's negative binomial half. A candidate is read by reweighting the default's run, or from
its own run where the weights fail. The default reading is what users get: today's sampler at `bart()`'s default
lengths.

Cell names: a real dataset as fitted (`roaches`), with a BART-fitted truth (`-bart`), widened to a spread of 2 or
3 (`-w2`, `-w3`), with a shape below 1 (`-r`), subsampled to 30 or 100 rows and widened (`-n30-3`), a designed
cell (`D1` to `D8`; D5 is null data), and a dataset's real outcomes over 5 folds (`-folds`).

## Part A, the prior alone

Percent of prior datasets. A1, spread above 4, at the worst design, at most 1 to pass. A2, a count the package
refuses, at the worst design, at most 1 to pass. A3, each real dataset no rarer than 1 in 100, at crimtab, at
least 1 to pass. The maintainer ruled that A3 is read at the middle of a dataset's four fits' range, 2.27 on
crimtab.

| prior | A1 | A2 | A3 at 2.27 | verdict |
|---|---|---|---|---|
| C0, chi(1.5, 2), the default | 8.41 | 0.77 | 8.90 | fails A1 |
| C2 chi(3, 1.25) | 2.00 | 0.06 | 3.18 | fails A1 |
| C2 chi(4, 1.03) | 0.99 | 0.04 | 1.71 | passes |
| C2 chi(5, 0.9) | 0.52 | 0 | 1.15 | passes |
| C2 chi(7, 0.75) | 0.17 | 0 | 0.57 | fails A3 |
| C2 chi(10, 0.62) | 0.04 | 0 | 0.26 | fails A3 |
| C3 chi(1.5, 3) | 4.62 | 0.43 | 5.02 | fails A1 |
| C4, k at least 0.6 | 0.89 | 0 | 2.14 | passes |
| C4, k at least 0.75 | 0.27 | 0 | 1.15 | passes |
| C4, k at least 1 | 0.02 | 0 | 0.33 | fails A3 |
| C5, half-normal, median 1.3 | 0.95 | 0 | 2.07 | passes |
| C5, half-normal, median 1.57 | 2.29 | 0.01 | 3.96 | fails A1 |
| C1, k = 1.5 | 0 | 0 | 0.19 | fails A3 |
| C1, k = 2 | 0 | 0 | 0 | fails A3 |
| C7, at each real design's own spread | 7.72 | 0 | 0.02 | fails A1 and A3 |

Five pass, none with margin, and go to the fits. The default fails, so step 4 of the rule never applies.

## What the first look decided

Look 1 read 540 replications, 485 on truths and 55 folds. Its information fraction is 0.379 and its boundary
3.46: a group-sequential statistic passes at z of at least 3.46 and fails at z of at most -3.46. Step 2 fails if
any part fails. A look decides when the recommendation is the same under every way its pending results could
resolve. All five fail, so only one outcome is left: the default. Step 3 has no pair to read, and step 5 finds no
candidate to check, so no guard runs and no block is funded. This look decides because no pending result can change
the outcome. It is not what the plan calls a final look, where pending results are read as failures and the
information fraction is set to 1.

Each candidate's gap less the default's (standard error) and status. The gap is the distance from the expected
coverage on either side, corrected for noise. The bar is 0.02.

| statistic | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 |
|---|---|---|---|---|---|
| F | -0.002 (0.001) pass | -0.003 (0.001) pass | -0.001 (0.001) pass | -0.002 (0.001) pass | +0.002 (0.002) pass |
| G | +0.007 (0.002) pass | +0.006 (0.003) pass | +0.009 (0.002) pass | +0.016 (0.002) pending | +0.010 (0.003) pending |
| R | -0.007 (0.003) pass | -0.010 (0.003) pass | -0.011 (0.002) pass | -0.005 (0.004) pass | -0.006 (0.015) pending |
| P | -0.000 (0.001) pass | -0.000 (0.001) pass | +0.000 (0.001) pass | -0.001 (0.001) pass | +0.000 (0.001) pass |
| widest F | -0.001 (0.002) pass | -0.001 (0.003) pass | -0.005 (0.003) pass | -0.004 (0.004) pass | -0.003 (0.003) pass |
| widest G | +0.006 (0.006) pending | +0.008 (0.004) pending | +0.022 (0.008) pending | +0.026 (0.006) pending | +0.019 (0.005) pending |
| widest R | -0.016 (0.010) pass | -0.016 (0.010) pass | -0.016 (0.011) pending | -0.014 (0.018) pending | -0.016 (0.010) pass |
| widest P | +0.000 (0.001) pass | +0.000 (0.001) pass | -0.002 (0.003) pass | -0.003 (0.003) pass | +0.001 (0.002) pass |

Interval score, candidate over default (standard error) and status. The bar is 1.05.

| statistic | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 |
|---|---|---|---|---|---|
| interval score F | 0.76 (0.022) pass | 0.73 (0.023) pass | 0.64 (0.027) pass | 0.61 (0.028) pass | 0.67 (0.026) pass |
| interval score G | 1.01 (0.029) pending | 1.04 (0.032) pending | 1.25 (0.043) fail | 1.39 (0.052) fail | 1.12 (0.032) pending |
| interval score R | 0.77 (0.018) pass | 0.74 (0.019) pass | 0.68 (0.026) pass | 0.65 (0.027) pass | 0.70 (0.025) pass |

The log score is C0's less the candidate's per row, the bar 0.01 at its 95 percent upper bound; all pass. Rank is
the mean of F's and R's gaps less C0's.

| row | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 |
|---|---|---|---|---|---|
| log score, C0 less candidate (95 upper) | +0.0030 (+0.0041) | +0.0041 (+0.0053) | -0.0003 (+0.0005) | +0.0004 (+0.0012) | +0.0008 (+0.0039) |
| rank | -0.0047 | -0.0065 | -0.0061 | -0.0035 | -0.0018 |

The width guard on each narrow cell: F's interval score ratio, candidate over default (95 percent interval by t)
and status. The bar is 1.10. C4 equals the default exactly on narrow cells, whose k stays above the floor.

| cell | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 |
|---|---|---|---|---|---|
| D5 | 1.21 (1.20, 1.23) fail | 1.26 (1.24, 1.28) fail | 1.00 (1.00, 1.00) pass | 1.00 (1.00, 1.00) pass | 0.41 (0.33, 0.48) pass |
| D8 | 1.03 (1.02, 1.03) pass | 1.03 (1.02, 1.04) pass | 1.00 (1.00, 1.00) pass | 1.00 (1.00, 1.00) pass | 0.99 (0.98, 1.00) pass |
| Insurance | 1.11 (1.09, 1.12) pending | 1.14 (1.12, 1.16) fail | 1.00 (1.00, 1.00) pass | 1.00 (1.00, 1.00) pass | 0.94 (0.89, 0.99) pass |
| Insurance-bart | 1.11 (1.09, 1.13) pending | 1.14 (1.11, 1.17) fail | 1.00 (1.00, 1.00) pass | 1.00 (1.00, 1.00) pass | 0.92 (0.89, 0.96) pass |
| Traffic | 1.16 (1.12, 1.19) fail | 1.20 (1.15, 1.24) fail | 1.00 (1.00, 1.00) pass | 1.00 (1.00, 1.00) pass | 1.03 (0.41, 1.65) pending |
| Traffic-bart | 0.92 (0.85, 0.99) pass | 0.90 (0.82, 0.98) pass | 1.00 (1.00, 1.00) pass | 1.00 (1.00, 1.00) pass | 1.41 (1.05, 1.78) pending |
| quakes | 1.04 (1.03, 1.06) pass | 1.06 (1.04, 1.08) pass | 1.00 (1.00, 1.00) pass | 1.00 (1.00, 1.00) pass | 0.98 (0.98, 0.99) pass |
| quakes-bart | 1.02 (1.01, 1.03) pass | 1.03 (1.03, 1.03) pass | 1.00 (1.00, 1.00) pass | 1.00 (1.00, 1.00) pass | 0.99 (0.98, 1.00) pass |

Cells worse by the bar, each candidate's gap less the default's on that cell:

| candidate | class | cell | difference | standard error | default's coverage | candidate's coverage |
|---|---|---|---|---|---|---|
| chi(4, 1.03) | G | epil-n30-3 | +0.219 | 0.050 | 0.733 | 0.517 |
| chi(5, 0.9) | G | epil-n30-3 | +0.252 | 0.051 | 0.733 | 0.483 |
| chi(5, 0.9) | G | epil-n100-3 | +0.101 | 0.027 | 0.783 | 0.683 |
| C4 at 0.6 | G | epil-n30-3 | +0.318 | 0.076 | 0.733 | 0.417 |
| C4 at 0.6 | G | epil-n100-3 | +0.199 | 0.065 | 0.783 | 0.583 |
| C4 at 0.75 | G | epil-n30-3 | +0.385 | 0.075 | 0.733 | 0.350 |
| C4 at 0.75 | G | epil-n100-3 | +0.250 | 0.067 | 0.783 | 0.533 |
| C4 at 0.75 | G | epil-w3 | +0.235 | 0.041 | 0.767 | 0.533 |
| C4 at 0.75 | G | D4-3 | +0.352 | 0.072 | 0.833 | 0.517 |
| C5 at 1.3 | G | epil-n30-3 | +0.303 | 0.065 | 0.733 | 0.433 |
| C5 at 1.3 | G | epil-w3 | +0.238 | 0.041 | 0.767 | 0.533 |
| C5 at 1.3 | F | Traffic-bart | +0.183 | 0.032 | 0.778 | 0.595 |

No shape-below-1 cell is worse by the bar for any candidate.

The flags on epil's three cells are a loss at the bottom of the range. Group readings on those cells that the default
covers and the candidate misses, of 150 in all:

| candidate | lost | true mean below 0.01 | truth below the interval | median true mean | median lower end, count scale |
|---|---|---|---|---|---|
| chi(4, 1.03) | 19 | 17 | 18 | 0.00019 | 0.0007 |
| chi(5, 0.9) | 25 | 20 | 22 | 0.00020 | 0.0013 |
| C4 at 0.6 | 42 | 35 | 40 | 0.00033 | 0.0032 |
| C4 at 0.75 | 49 | 36 | 47 | 0.00038 | 0.0089 |
| C5 at 1.3 | 35 | 31 | 33 | 0.00031 | 0.0015 |

On those cells groups at or above 0.01 are covered 0.82 by the default, and 0.73 to 0.82 by the candidates (C4 at
0.6, 0.80). Groups below 0.01 are covered 0.66 by the default, and 0 to 0.35 by the candidates (C4 at 0.6, 0.02).

How the five fail, and how firmly:

| candidate | inside the spending | outside it, fixed-level t tests |
|---|---|---|
| chi(4, 1.03) | pending | width guard on D5 and Traffic; one cell |
| chi(5, 0.9) | pending | width guard on D5, Insurance, Insurance-bart and Traffic; two cells |
| C4 at 0.6 | fails: G's interval score, z -4.57 | two cells |
| C4 at 0.75 | fails: G's interval score, z -6.52 | four cells |
| C5 at 1.3 | pending | three cells |

The cell rule and the width guard are outside the spending and repeat at each look, as the pre-registration
discloses. It also discloses that the cell rule reads high near the bar: a true difference of 0.021 flagged 14 of
68 cells in a constructed case, and 0.015 none. That does not touch these flags: the smallest flagged
difference is 0.101, five times the bar. The nearest width failure is Traffic under chi(4, 1.03), at 1.16 with a
lower bound of 1.12 against 1.10.

## Readings by class and band

90 percent coverage on the long reading, cells equal within a band, against the expected coverage (0.9 N / (N + 1)
at the run's effective draws; 0.9 for P).

F:

| band | C0 | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 | expected |
|---|---|---|---|---|---|---|---|
| narrow | 0.947 | 0.952 | 0.954 | 0.947 | 0.947 | 0.909 | 0.899 |
| middle | 0.926 | 0.926 | 0.926 | 0.926 | 0.926 | 0.927 | 0.898 |
| wide | 0.931 | 0.931 | 0.931 | 0.933 | 0.933 | 0.933 | 0.898 |
| small | 0.888 | 0.893 | 0.894 | 0.897 | 0.900 | 0.895 | 0.899 |

G:

| band | C0 | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 | expected |
|---|---|---|---|---|---|---|---|
| narrow | 0.931 | 0.927 | 0.929 | 0.931 | 0.931 | 0.921 | 0.899 |
| middle | 0.848 | 0.848 | 0.855 | 0.848 | 0.852 | 0.848 | 0.899 |
| wide | 0.891 | 0.885 | 0.885 | 0.880 | 0.873 | 0.879 | 0.899 |
| small | 0.838 | 0.822 | 0.818 | 0.810 | 0.782 | 0.811 | 0.899 |

R:

| band | C0 | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 | expected |
|---|---|---|---|---|---|---|---|
| narrow | 0.897 | 0.907 | 0.906 | 0.897 | 0.897 | 0.871 | 0.899 |
| middle | 0.892 | 0.892 | 0.892 | 0.892 | 0.892 | 0.893 | 0.898 |
| wide | 0.961 | 0.954 | 0.948 | 0.954 | 0.947 | 0.954 | 0.898 |
| small | 0.971 | 0.960 | 0.955 | 0.934 | 0.913 | 0.942 | 0.899 |

P:

| band | C0 | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 | expected |
|---|---|---|---|---|---|---|---|
| narrow | 0.921 | 0.921 | 0.920 | 0.921 | 0.921 | 0.916 | 0.900 |
| middle | 0.888 | 0.888 | 0.888 | 0.888 | 0.888 | 0.888 | 0.900 |
| wide | 0.897 | 0.897 | 0.896 | 0.898 | 0.899 | 0.897 | 0.900 |
| small | 0.900 | 0.892 | 0.892 | 0.897 | 0.895 | 0.897 | 0.900 |

The 50, 80 and 95 percent coverages, read beside and not by the rule; each entry gives the three in that order.
P is read at 90 percent only.

F:

| band | C0 | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 |
|---|---|---|---|---|---|---|
| narrow | 0.58 / 0.88 / 0.97 | 0.61 / 0.89 / 0.98 | 0.61 / 0.89 / 0.98 | 0.58 / 0.88 / 0.97 | 0.58 / 0.88 / 0.97 | 0.54 / 0.83 / 0.95 |
| middle | 0.55 / 0.84 / 0.96 | 0.55 / 0.84 / 0.96 | 0.55 / 0.84 / 0.96 | 0.55 / 0.84 / 0.96 | 0.55 / 0.84 / 0.96 | 0.55 / 0.85 / 0.96 |
| wide | 0.55 / 0.84 / 0.97 | 0.55 / 0.84 / 0.97 | 0.55 / 0.84 / 0.97 | 0.55 / 0.85 / 0.97 | 0.56 / 0.85 / 0.97 | 0.55 / 0.85 / 0.97 |
| small | 0.51 / 0.78 / 0.94 | 0.52 / 0.80 / 0.94 | 0.52 / 0.80 / 0.94 | 0.52 / 0.81 / 0.96 | 0.51 / 0.80 / 0.95 | 0.52 / 0.80 / 0.95 |

G:

| band | C0 | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 |
|---|---|---|---|---|---|---|
| narrow | 0.55 / 0.85 / 0.96 | 0.59 / 0.87 / 0.95 | 0.56 / 0.87 / 0.96 | 0.55 / 0.85 / 0.96 | 0.55 / 0.85 / 0.96 | 0.53 / 0.81 / 0.94 |
| middle | 0.45 / 0.73 / 0.91 | 0.44 / 0.73 / 0.91 | 0.45 / 0.73 / 0.91 | 0.45 / 0.73 / 0.91 | 0.44 / 0.73 / 0.91 | 0.44 / 0.72 / 0.91 |
| wide | 0.46 / 0.76 / 0.94 | 0.47 / 0.76 / 0.94 | 0.47 / 0.75 / 0.94 | 0.47 / 0.75 / 0.92 | 0.46 / 0.75 / 0.92 | 0.46 / 0.75 / 0.93 |
| small | 0.44 / 0.73 / 0.91 | 0.42 / 0.72 / 0.88 | 0.42 / 0.71 / 0.88 | 0.40 / 0.69 / 0.87 | 0.38 / 0.67 / 0.86 | 0.41 / 0.71 / 0.88 |

R:

| band | C0 | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 |
|---|---|---|---|---|---|---|
| narrow | 0.45 / 0.79 / 0.94 | 0.51 / 0.83 / 0.97 | 0.54 / 0.83 / 0.97 | 0.45 / 0.79 / 0.94 | 0.45 / 0.79 / 0.94 | 0.50 / 0.77 / 0.91 |
| middle | 0.48 / 0.77 / 0.97 | 0.48 / 0.78 / 0.97 | 0.48 / 0.79 / 0.97 | 0.48 / 0.78 / 0.97 | 0.48 / 0.77 / 0.97 | 0.48 / 0.78 / 0.97 |
| wide | 0.60 / 0.90 / 0.99 | 0.61 / 0.90 / 0.99 | 0.60 / 0.90 / 0.99 | 0.59 / 0.89 / 0.99 | 0.58 / 0.87 / 0.98 | 0.59 / 0.89 / 0.99 |
| small | 0.61 / 0.89 / 0.99 | 0.60 / 0.88 / 0.99 | 0.60 / 0.88 / 0.99 | 0.55 / 0.84 / 0.98 | 0.54 / 0.82 / 0.96 | 0.57 / 0.85 / 0.98 |

The shape r is reported and not read. The median of the default's runs' median r, by group: narrow 17.5, middle 3,
wide 1, small 8, the shape-below-1 cells 1 (the grid's floor), the folds 12.

## Rows whose true mean is below 0.01

Scored apart and not read by the rule. Their coverage is set by how far each prior lets the log mean fall. Long
reading, F's rows and then P's. Only the middle, wide and small bands hold such rows.

F's rows:

| band | C0 | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 |
|---|---|---|---|---|---|---|
| middle | 0.654 | 0.647 | 0.653 | 0.652 | 0.639 | 0.652 |
| wide | 0.849 | 0.830 | 0.821 | 0.768 | 0.723 | 0.804 |
| small | 0.886 | 0.817 | 0.786 | 0.684 | 0.570 | 0.753 |

P's rows:

| band | C0 | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 |
|---|---|---|---|---|---|---|
| middle | 0.892 | 0.892 | 0.892 | 0.892 | 0.892 | 0.892 |
| wide | 0.905 | 0.908 | 0.908 | 0.908 | 0.909 | 0.908 |
| small | 0.930 | 0.936 | 0.936 | 0.937 | 0.945 | 0.932 |

## The default reading

Guard readings were not reached, so the default reading exists for the default only. The default against itself,
long and default length, with the mean width on the log scale:

| class | band | long cover | long gap | default cover | default gap | long width | default width |
|---|---|---|---|---|---|---|---|
| F | narrow | 0.947 | 0.078 | 0.943 | 0.078 | 0.43 | 0.43 |
| F | middle | 0.926 | 0.038 | 0.918 | 0.032 | 1.08 | 1.04 |
| F | wide | 0.931 | 0.038 | 0.921 | 0.031 | 1.62 | 1.45 |
| F | small | 0.888 | 0.036 | 0.894 | 0.032 | 5.65 | 3.07 |
| G | narrow | 0.931 | 0.037 | 0.931 | 0.037 | 0.17 | 0.17 |
| G | middle | 0.848 | 0.083 | 0.840 | 0.086 | 0.40 | 0.39 |
| G | wide | 0.891 | 0.039 | 0.871 | 0.046 | 0.66 | 0.58 |
| G | small | 0.838 | 0.072 | 0.809 | 0.096 | 2.39 | 1.53 |
| R | narrow | 0.897 | 0.068 | 0.894 | 0.069 | 0.35 | 0.34 |
| R | middle | 0.892 | 0.088 | 0.868 | 0.094 | 0.75 | 0.72 |
| R | wide | 0.961 | 0.086 | 0.935 | 0.062 | 1.85 | 1.50 |
| R | small | 0.971 | 0.067 | 0.941 | 0.041 | 6.58 | 3.59 |
| P | narrow | 0.921 | 0.023 | 0.921 | 0.023 | - | - |
| P | middle | 0.888 | 0.018 | 0.888 | 0.019 | - | - |
| P | wide | 0.897 | 0.011 | 0.897 | 0.010 | - | - |
| P | small | 0.900 | 0.034 | 0.897 | 0.034 | - | - |
| Flow | narrow | - | - | - | - | - | - |
| Flow | middle | 0.654 | 0.267 | 0.617 | 0.290 | 1.99 | 1.85 |
| Flow | wide | 0.848 | 0.142 | 0.725 | 0.236 | 4.49 | 3.61 |
| Flow | small | 0.886 | 0.089 | 0.726 | 0.200 | 12.92 | 6.16 |
| Plow | narrow | - | - | - | - | - | - |
| Plow | middle | 0.892 | 0.056 | 0.892 | 0.056 | - | - |
| Plow | wide | 0.908 | 0.021 | 0.911 | 0.016 | - | - |
| Plow | small | 0.930 | 0.029 | 0.931 | 0.029 | - | - |

At default length k has not mixed, and the intervals on wide small-sample data are about half as wide as the long
reading's (F 3.07 against 5.65). Coverage on F and P holds; on G and R it falls by 0.02 to 0.03 in the wide and
small bands, and on R in the middle band too, and on the rows below 0.01 it falls from 0.85 to 0.73 on wide data and from 0.89 to 0.73 on small.

If the collapsed move later lands for negative binomial, the default reading changes and the long reading does not.

## What keeping the default costs, and what a change would cost

**Keeping it.** The default keeps a prior whose top tail spends its mass where no count data have been seen: 4.0 to
8.4 percent of prior datasets past a spread of 4, against the 1 percent of Part A. A change of default was expected
to buy a coherent prior predictive, not better intervals, and that stays unbought.

It also keeps some poor intervals on the widest small samples. On D4-3 (n = 60, spread 3) F covers 0.54 with a gap of
0.356, and a mean width of 46 on the log scale; the candidates cover 0.63 to 0.77, their gaps smaller by 0.088 to
0.227. F's interval score there is 64.3 under the default against 6.7 to 24.0. That one cell is 68 to 76 percent of
the fall in F's interval score: 0.69 of 0.97 for C4 at 0.6. Over the small band F's score is 7.18 under the default
and 3.22 to 4.57 under the candidates, and over the wide band 1.90 against 1.64 to 1.85. RMSE of the log mean moves
with it (wide band 0.494 against 0.398 to 0.476).

The default also keeps weak coverage of group means on epil's widened cells: 0.73 to 0.78 on the long reading against
an expected 0.90, and 0.48 to 0.63 at default length, which is what users get. The candidates' long readings on those
cells cover 0.35 to 0.75, so at default length the default sits about where they do on the long reading. Their default
readings were not run, so the two cannot be compared like with like.

**Changing it.** Every candidate gives something up:

- the chi priors, wider intervals on narrow data, where the default does well (interval score 1.21 and 1.26 times the
  default's on null data);
- every candidate, coverage on rows and groups whose true mean is near zero. On rows below 0.01 F's coverage falls
  from 0.89 to between 0.57 and 0.82 on small data, and from 0.85 to between 0.72 and 0.83 on wide;
- C5, coverage on Traffic-bart (0.60 against 0.78); C4 at 0.75, coverage on D4-3's group means (0.52 against 0.83).

## The long reading's check, run lengths and what stayed unsettled

**The check.** Where both builds settled on stage 2's seven cells, k's 5th, 50th and 95th percentiles agree between
the scratch build and today's sampler (21 z's, chi-squared 28.0 on 21, p 0.14, the largest |z| 2.55), and the
rows' interval endpoints agree within two combined standard errors in 92 to 97 percent of rows. On the run machine
the three libraries passed the identity check, 8 of 8 pairs.

**Lengths and own runs.** Sweeps a chain of the default's run; own runs a replication, of five fitted candidates.

| group | replications | median sweeps | largest | own runs a replication | runs unsettled | of runs |
|---|---|---|---|---|---|---|
| narrow | 40 | 2,500 | 5,000 | 2.77 | 0 | 151 |
| middle | 70 | 2,500 | 15,000 | 0.27 | 0 | 89 |
| wide | 145 | 7,500 | 100,000 | 1.77 | 18 | 401 |
| small | 200 | 5,000 | 77,500 | 2.33 | 0 | 667 |
| shape below 1 | 30 | 2,500 | 20,000 | 1.17 | 0 | 65 |
| folds | 55 | 5,000 | 100,000 | 2.00 | 30 | 165 |

**Unsettled.** 48 of 1,538 runs reached the cap of 100,000 sweeps short of a target: effective draws of the rows
53 to 99, of k 74 and up. All are on `mite-bart` (a run on each of its 5 replications; the default's own on 3)
and on mite's 5 folds. Those two cells pass the rule's 20 percent, so the rule reports them apart from their band's
averages; one wide cell of 26 does not reach the third of a band that would stop the band's reading. The harness
kept them in the band averages; see [Disclosures](#disclosures). Taken out, every candidate's verdict and the
outcome are unchanged.

**Mixing, for the collapsed move's plan.** Total sweeps per effective draw of k, at each run's stop. The same datasets
under today's sampler and under the scratch build, on the real outcomes. Today's figures and the first column for
the move are the probes' long runs; the study's are the median of 5 folds, each four fifths of the rows; mite's are at
the cap.

| real outcome | today's sampler, probes | the move, probes | the move, this study's folds | today over this study |
|---|---|---|---|---|
| Traffic | 620 | 37 | 24 | 26 times |
| roaches | 1,200 | 236 | 146 | 8 times |
| solder | 3,720 | 1,430 | 1,403 | 2.7 times |
| mite | 8,000 | 1,600 | 2,010 | 4 times |

The move costs 17 to 24 percent more a sweep, so the gain in effective draws a second is about a fifth less than
these ratios. Over all 11 datasets' folds the median is 18 (Insurance) to 2,010 (mite). The gain on k does not reach
the rows on the slowest datasets: the probes' slowest estimand went from 1,910 to 2,390 on solder and from 6,110 to
4,020 on mite, and on mite's folds the rows, not k, set the length, at 60 effective draws where k has 199.

Counts simulated from truths mix faster than real outcomes do, and do not compare with the table: the median is 28 on
narrow, 89 on middle, 239 on wide and 157 on small. The rule stops on effective draws and not on R-hat: on the default's
runs k's split R-hat at the stop has median 1.035, 90th percentile 1.09 and largest 1.55 (1.68 over all runs), and the
mean count's largest is 1.09.

## The shape's prior

C6 is a separate decision. The look scored it on only 133 of 540 replications, and the reason is a harness defect,
not the weights. The weights take at most 13 values, one for each point on the shape's grid, and the Pareto tail fit
usually returns no value on so few distinct values. The job counted that as a failure. So k-hat is missing on 359 of 540
replications and above 0.7 on 15 more, while the weights keep at least 47 percent of the draws on all 540 (median 89
percent). The decision is not touched, since C6 is not among the five fitted candidates.

Read on every replication from the 500 draws each run stores, with the harness's own reading and statistics and no
fit:

| statistic | C6 less the default, or C6 over the default | status at 3.46 |
|---|---|---|
| F | -0.0022 (0.0008) | pass |
| G | -0.0027 (0.0026) | pass |
| R | +0.0055 (0.0025) | pass |
| P | -0.0006 (0.0009) | pass |
| widest F | -0.0034 (0.0024) | pass |
| widest G | +0.0172 (0.0116) | pending |
| widest R | +0.0004 (0.0001) | pass |
| widest P | +0.0001 (0.0007) | pass |
| interval score F | 0.999 (0.007) | pass |
| interval score G | 0.995 (0.005) | pass |
| interval score R | 1.000 (0.004) | pass |

No cell and no shape-below-1 cell is worse by the bar, the width guard passes on all 8 narrow cells (0.97 to 1.00),
and the log score is 0.0017 nats a row better than the default's. Thinning costs little: a replication's coverage from
500 draws differs from its full reading by 0.008 on average, for both the default and C6. C6 leaves the spread's prior
as it is, so it does nothing for Part A. Step 2 reads pending on widest G alone, mostly one cell, quakes-w3 (0.85 to
0.75 over 5 replications). On this reading it neither loses nor gains much.

## Cost

Look 1 in box core-hours, against the pre-registered sizing at the calibrated speed. In the sizing table's units,
in which the stop is set, it comes to 94.5 of 110. By component, in the harness's seconds, which the stop reads:

| component | sized | measured |
|---|---|---|
| the default's long runs | 12.67 | 21.18 |
| own runs | 16.41 | 42.97 |
| readings from all draws | 1.83 | 1.44 |
| the default's default readings | 1.55 | 0.96 |
| all | 32.46 | 66.55 |

By group, in the workers' seconds, which add R's start and the saving of each job's output:

| group | sized | measured |
|---|---|---|
| narrow | 3.31 | 3.73 |
| middle | 1.90 | 3.39 |
| wide | 18.45 | 32.87 |
| small | 5.28 | 7.88 |
| shape below 1 | 0.92 | 1.97 |
| folds | 2.54 | 16.92 |
| all | 32.39 | 66.76 |

Almost all of the overrun is real price. The workers' processes took 238,942 CPU seconds against 240,342 s of job
time, so waiting for a core cost at most 0.6 percent. Other guests on the same machine added about 0.15 to its load for the first seven
hours and about 0.44 for the next nine, and can have slowed the workers through shared cores: at most about 2.5
core-hours at a 50 percent penalty on the core shared. Fixed-length work shows nothing of the kind: the default's
readings ran under their price, and the narrow band's 1.13 times its price is its own-run count, 3.77 against 3.37
runs a replication. The calibration's seconds, on the three cells it timed, are predicted from the main runs' own
within 1.00, 1.13 and 0.91, once its shorter burn-in is allowed for.

The sizing missed in three ways. Own runs fall where runs are longest: five own runs a replication on mite-w3,
epil-w3, solder-w3 and mite's folds, 4.4 on mite-bart, against band rates of 1.25 to 2.37 applied to average cost.
The real-outcome folds were priced from cells where runs were short, and mite's five hit the cap with five own runs
each, 10.4 core-hours against 2.5 for all 55 folds. Stage 2's cells did not include mite-bart or the folds, and it
reported that no slow-tree dataset neared the cap; 48 runs at the cap took 15.6 core-hours.

Each candidate's own runs, which are 0.56 to 0.92 of the default's length on the same replication:

| candidate | own runs | core-hours | seconds a run |
|---|---|---|---|
| chi(4, 1.03) | 206 | 11.07 | 193 |
| chi(5, 0.9) | 264 | 12.49 | 170 |
| C4 at 0.6 | 138 | 5.03 | 131 |
| C4 at 0.75 | 198 | 7.01 | 127 |
| C5 at 1.3 | 192 | 7.38 | 138 |
| the default's long run, for comparison | 540 | 21.18 | 141 |

## Added after the first look: two chi priors with a raised scale

This section was added after the first look, at the maintainer's question whether other settings of the chi prior had
been tried. Its two candidates were not pre-registered, and they were chosen after a rough look at the same stored runs
it reads. The look's outcome as registered is unchanged: the rule keeps the default, and the five registered candidates
fail. What follows is read under the same rule, and is evidence about those two priors, not a new outcome of the study.

**The result.** Both priors pass the prior check, without margin. On the fits neither fails any part of the rule, and
both lose far less near zero than the five did. Two readings of group means are still pending for each at the first
look's boundary, so at that boundary the rule cannot say. This reading has used up what the stop allowed, and under the rule a
look that can fund no more runs is final. Read as final, both pass and tie, the order is chi(2.75, 2) then chi(3, 1.75),
and chi(2.75, 2) also passes the default-length check. So read, the rule would choose chi(2.75, 2) over the default. That
finality exists only because this added reading spent the room. No fresh data had tested either prior at this point;
the next section reports the test of chi(2.75, 2).

### Why these two

Every chi prior tried but one held the default's median k at 1.9 and raised the degrees of freedom. That takes mass away
near zero, which was the aim, and also from the top, which was not. The tried priors failed mainly for the second
reason: on null data k's posterior lies above 4, where they put little mass. Raising the scale takes mass from near zero
and adds to the top. It also moves the middle: the median of k goes from 1.9 to about 2.8, so this changes the bulk of the
prior and not only its tail.

| prior | median k | median of 3/k | P(k < 0.5) | P(k < 1) | P(k > 4) |
|---|---|---|---|---|---|
| chi(1.5, 2), the default | 1.91 | 1.57 | 0.080 | 0.217 | 0.085 |
| chi(4, 1.03) | 1.89 | 1.59 | 0.006 | 0.082 | 0.005 |
| chi(5, 0.9) | 1.88 | 1.60 | 0.003 | 0.058 | 0.001 |
| chi(3, 1.75) | 2.69 | 1.11 | 0.006 | 0.045 | 0.156 |
| chi(2.75, 2) | 2.91 | 1.03 | 0.007 | 0.044 | 0.227 |

### Part A

Percent of prior datasets, standard errors in parentheses, with the bars of [Part A, the prior alone](#part-a-the-prior-alone): A1 and A2
at most 1, A3 at least 1, A3 read at crimtab's 2.27.

| prior | A1 | A2 | A3 at 2.27 | verdict |
|---|---|---|---|---|
| chi(3, 1.75) | 0.77 (0.013) | 0.04 (0.016) | 1.00 (0.070), and 1.12 (0.016) on 21 times the forests | passes |
| chi(2.75, 2) | 0.83 (0.008) | 0.04 (0.013) | 1.14 (0.075), and 1.17 (0.017) on 21 times the forests | passes |

Neither has margin. A3 at the planned number of forests sat on the bar for chi(3, 1.75), which is why it was rerun on
more. The same code on new forests gives 8.41, 0.82 and 9.21 for the default and 1.03, 0.04 and 1.73 for chi(4, 1.03),
close to the table above. Read on the forests stored from the first look, which reproduce that table's A1, the two
priors' A1 are 0.75 and 0.82.

### How they were read

The registered candidates were read by reweighting the default's long runs where the weights hold, and otherwise by
their own runs. The same is done here, under the same bounds and boundaries, with one difference. The draws that the
first look reweighted were not kept, only 500 from each run. So on a replication with no own run, a new prior's reading
is the default's registered reading plus the difference between the two priors on those 500 draws.

The first expectation, that a chi prior's weights would hold nearly everywhere, was wrong. The weights' tail check and
share of draws never come near failing (k-hat at most 0.31, never missing). But the rule also asks for 100 effective
draws, and the default's runs stop at a median of about 124, so little room is left: own runs are required on 116 and 93 of the 540
replications. 181 of these 209 ran, along with 50 more on narrow cells as a check. The narrow checks agree with the
reweighted readings.

As a control, the same composition was applied to the five registered candidates. It gives all their step 2 verdicts,
flagged cells and width statuses as registered. The largest difference in an estimate is 0.007, and three statistics
near the boundary move between pass and pending. A pending near the boundary is therefore not a firm reading here.

### Step 2, beside the five

Gap less the default's (standard error) and status at the boundary of 3.46. The bar is 0.02.

| statistic | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 | chi(3, 1.75) | chi(2.75, 2) |
|---|---|---|---|---|---|---|---|
| F | -0.002 (0.001) pass | -0.003 (0.001) pass | -0.001 (0.001) pass | -0.002 (0.001) pass | +0.002 (0.002) pass | -0.001 (0.001) pass | +0.000 (0.001) pass |
| G | +0.007 (0.002) pass | +0.006 (0.003) pass | +0.009 (0.002) pass | +0.016 (0.002) pending | +0.010 (0.003) pending | +0.003 (0.002) pass | +0.004 (0.002) pass |
| R | -0.007 (0.003) pass | -0.010 (0.003) pass | -0.011 (0.002) pass | -0.005 (0.004) pass | -0.006 (0.015) pending | -0.002 (0.002) pass | -0.006 (0.003) pass |
| P | -0.000 (0.001) pass | -0.000 (0.001) pass | +0.000 (0.001) pass | -0.001 (0.001) pass | +0.000 (0.001) pass | -0.001 (0.001) pass | -0.000 (0.001) pass |
| widest F | -0.001 (0.002) pass | -0.001 (0.003) pass | -0.005 (0.003) pass | -0.004 (0.004) pass | -0.003 (0.003) pass | -0.001 (0.001) pass | -0.001 (0.001) pass |
| widest G | +0.006 (0.006) pending | +0.008 (0.004) pending | +0.022 (0.008) pending | +0.026 (0.006) pending | +0.019 (0.005) pending | +0.005 (0.004) pending | +0.008 (0.004) pending |
| widest R | -0.016 (0.010) pass | -0.016 (0.010) pass | -0.016 (0.011) pending | -0.014 (0.018) pending | -0.016 (0.010) pass | -0.000 (0.000) pass | -0.017 (0.010) pass |
| widest P | +0.000 (0.001) pass | +0.000 (0.001) pass | -0.002 (0.003) pass | -0.003 (0.003) pass | +0.001 (0.002) pass | +0.000 (0.001) pass | +0.000 (0.001) pass |

Interval score, candidate over default (standard error) and status. The bar is 1.05.

| statistic | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 | chi(3, 1.75) | chi(2.75, 2) |
|---|---|---|---|---|---|---|---|
| interval score F | 0.76 (0.022) pass | 0.73 (0.023) pass | 0.64 (0.027) pass | 0.61 (0.028) pass | 0.67 (0.026) pass | 0.80 (0.020) pass | 0.82 (0.019) pass |
| interval score G | 1.01 (0.029) pending | 1.04 (0.032) pending | 1.25 (0.043) fail | 1.39 (0.052) fail | 1.12 (0.032) pending | 0.99 (0.023) pending | 0.99 (0.022) pending |
| interval score R | 0.77 (0.018) pass | 0.74 (0.019) pass | 0.68 (0.026) pass | 0.65 (0.027) pass | 0.70 (0.025) pass | 0.81 (0.017) pass | 0.83 (0.015) pass |
| same, without D4-3: F | 0.92 | 0.90 | 0.85 | 0.82 | 0.87 | 0.94 | 0.94 |
| same, without D4-3: G | 1.10 | 1.15 | 1.40 | 1.54 | 1.25 | 1.06 | 1.06 |
| same, without D4-3: R | 0.89 | 0.87 | 0.82 | 0.79 | 0.85 | 0.91 | 0.92 |

The rows without D4-3 are read beside the rule. G's ratio of 0.99 is under its bar with that one design and just over it
without. By cell, the median G ratio is 1.00 for both priors. Above 1.05 are 5 and 4 of 68 cells, led by epil-n30-3 (1.30,
1.25), epil-n100-3 (1.27, 1.27) and epil-w3 (1.11, 1.12).

Log score, C0's less the candidate's per row (95 percent upper bound; the bar is 0.01), and rank, the mean of F's and R's
gaps less C0's:

| row | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 | chi(3, 1.75) | chi(2.75, 2) |
|---|---|---|---|---|---|---|---|
| log score | +0.0030 (+0.0041) | +0.0041 (+0.0053) | -0.0003 (+0.0005) | +0.0004 (+0.0012) | +0.0008 (+0.0039) | +0.0005 (+0.0009) | +0.0002 (+0.0007) |
| rank | -0.0047 | -0.0065 | -0.0061 | -0.0035 | -0.0018 | -0.0017 | -0.0030 |

The width guard, F's interval score ratio on each narrow cell (95 percent interval). The bar is 1.10. This is where the
tried chi priors failed.

| cell | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 and 0.75 | C5 at 1.3 | chi(3, 1.75) | chi(2.75, 2) |
|---|---|---|---|---|---|---|
| D5 | 1.21 (1.20, 1.23) fail | 1.26 (1.24, 1.28) fail | 1.00 pass | 0.41 (0.33, 0.48) pass | 1.02 (1.01, 1.03) pass | 0.98 (0.97, 0.98) pass |
| D8 | 1.03 (1.02, 1.03) pass | 1.03 (1.02, 1.04) pass | 1.00 pass | 0.99 (0.98, 1.00) pass | 1.00 (1.00, 1.00) pass | 0.99 (0.99, 1.00) pass |
| Insurance | 1.11 (1.09, 1.12) pending | 1.14 (1.12, 1.16) fail | 1.00 pass | 0.94 (0.89, 0.99) pass | 1.01 (1.00, 1.01) pass | 0.98 (0.97, 1.00) pass |
| Insurance-bart | 1.11 (1.09, 1.13) pending | 1.14 (1.11, 1.17) fail | 1.00 pass | 0.92 (0.89, 0.96) pass | 1.01 (1.00, 1.01) pass | 0.98 (0.97, 1.00) pass |
| Traffic | 1.16 (1.12, 1.19) fail | 1.20 (1.15, 1.24) fail | 1.00 pass | 1.03 (0.41, 1.65) pending | 1.02 (1.01, 1.03) pass | 0.99 (0.97, 1.01) pass |
| Traffic-bart | 0.92 (0.85, 0.99) pass | 0.90 (0.82, 0.98) pass | 1.00 pass | 1.41 (1.05, 1.78) pending | 1.00 (0.99, 1.01) pass | 1.03 (1.00, 1.05) pass |
| quakes | 1.04 (1.03, 1.06) pass | 1.06 (1.04, 1.08) pass | 1.00 pass | 0.98 (0.98, 0.99) pass | 1.00 (1.00, 1.00) pass | 0.99 (0.99, 1.00) pass |
| quakes-bart | 1.02 (1.01, 1.03) pass | 1.03 (1.03, 1.03) pass | 1.00 pass | 0.99 (0.98, 1.00) pass | 1.00 (1.00, 1.00) pass | 0.99 (0.99, 1.00) pass |

No cell is worse by the bar for either prior, and no shape-below-1 cell. The five had one to four flagged cells each.
That says little: at 5 replications a cell, a loss of the same size in each replication it touches is flagged only if it
touches at least 4 of the 5. The largest differences are still large, and they escape the flag for that reason: G on epil-n30-3,
+0.086 (0.045) and +0.052 (0.043), and on epil-w3, +0.070 (0.041) and +0.105 (0.041). The nearest to a flag is F on
Traffic-bart under chi(2.75, 2), +0.026 (0.005), whose lower bound is 0.012 against 0.02.

How the seven stand:

| candidate | inside the spending, at 3.46 | outside it, fixed-level t tests |
|---|---|---|
| chi(4, 1.03) | pending | fails: width guard on D5 and Traffic; one cell |
| chi(5, 0.9) | pending | fails: width guard on D5, Insurance, Insurance-bart and Traffic; two cells |
| C4 at 0.6 | fails: G's interval score, z -4.57 | fails: two cells |
| C4 at 0.75 | fails: G's interval score, z -6.52 | fails: four cells |
| C5 at 1.3 | pending | fails: three cells |
| chi(3, 1.75) | pending: widest G z 3.33, G's interval score z 2.65 | passes |
| chi(2.75, 2) | pending: widest G z 2.68, G's interval score z 2.73 | passes |

### What the rule says

At the first look's boundary of 3.46 nothing decides. The recommendation could still be the default, chi(2.75, 2) or
chi(3, 1.75), and the two priors' lead over each other is pending as well.

The next step would be a block of 5 more replications. There is no room for one, so the look is final: the boundary
becomes 1.96, pending reads as a failure, and a pair pending reads as tied. Then every pending part above passes. The two
are tied on rank, with a lead of 0.0013 (0.0009) against 0.005. They need the same code change, and the tie goes to the
smaller interval score on narrow data, 0.992 against 1.007, so chi(2.75, 2) comes first.

### Coverage by band

90 percent coverage on the long reading, cells equal within a band, as in
[Readings by class and band](#readings-by-class-and-band).

| class | band | C0 | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 | chi(3, 1.75) | chi(2.75, 2) |
|---|---|---|---|---|---|---|---|---|---|
| F | narrow | 0.947 | 0.952 | 0.954 | 0.947 | 0.947 | 0.909 | 0.947 | 0.944 |
| F | middle | 0.926 | 0.926 | 0.926 | 0.926 | 0.926 | 0.927 | 0.924 | 0.923 |
| F | wide | 0.931 | 0.931 | 0.931 | 0.933 | 0.933 | 0.933 | 0.930 | 0.931 |
| F | small | 0.888 | 0.893 | 0.894 | 0.897 | 0.900 | 0.895 | 0.891 | 0.891 |
| G | narrow | 0.931 | 0.927 | 0.929 | 0.931 | 0.931 | 0.921 | 0.931 | 0.931 |
| G | middle | 0.848 | 0.848 | 0.855 | 0.848 | 0.852 | 0.848 | 0.845 | 0.843 |
| G | wide | 0.891 | 0.885 | 0.885 | 0.880 | 0.873 | 0.879 | 0.889 | 0.887 |
| G | small | 0.838 | 0.822 | 0.818 | 0.810 | 0.782 | 0.811 | 0.833 | 0.832 |
| R | narrow | 0.897 | 0.907 | 0.906 | 0.897 | 0.897 | 0.871 | 0.898 | 0.896 |
| R | middle | 0.892 | 0.892 | 0.892 | 0.892 | 0.892 | 0.893 | 0.882 | 0.881 |
| R | wide | 0.961 | 0.954 | 0.948 | 0.954 | 0.947 | 0.954 | 0.955 | 0.949 |
| R | small | 0.971 | 0.960 | 0.955 | 0.934 | 0.913 | 0.942 | 0.966 | 0.959 |
| P | narrow | 0.921 | 0.921 | 0.920 | 0.921 | 0.921 | 0.916 | 0.921 | 0.921 |
| P | middle | 0.888 | 0.888 | 0.888 | 0.888 | 0.888 | 0.888 | 0.888 | 0.888 |
| P | wide | 0.897 | 0.897 | 0.896 | 0.898 | 0.899 | 0.897 | 0.897 | 0.897 |
| P | small | 0.900 | 0.892 | 0.892 | 0.897 | 0.895 | 0.897 | 0.894 | 0.891 |

The 50, 80 and 95 percent coverages, read beside and not by the rule, stay within 0.03 of the default's in every class
and band. The largest move is R's 80 percent on narrow data, 0.79 to 0.76 under chi(3, 1.75).

### Rows and groups whose true mean is below 0.01

Scored apart and not read by the rule, for rows. This is where a prior that thins k near zero pays, and where the two
priors pay less than any of the five.

F's rows, long reading:

| band | C0 | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 | chi(3, 1.75) | chi(2.75, 2) |
|---|---|---|---|---|---|---|---|---|
| middle | 0.654 | 0.647 | 0.653 | 0.652 | 0.639 | 0.652 | 0.651 | 0.653 |
| wide | 0.849 | 0.830 | 0.821 | 0.768 | 0.723 | 0.804 | 0.836 | 0.839 |
| small | 0.886 | 0.817 | 0.786 | 0.684 | 0.570 | 0.753 | 0.845 | 0.854 |

P's rows, long reading:

| band | C0 | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 | chi(3, 1.75) | chi(2.75, 2) |
|---|---|---|---|---|---|---|---|---|
| middle | 0.892 | 0.892 | 0.892 | 0.892 | 0.892 | 0.892 | 0.892 | 0.892 |
| wide | 0.905 | 0.908 | 0.908 | 0.908 | 0.909 | 0.908 | 0.905 | 0.905 |
| small | 0.930 | 0.936 | 0.936 | 0.937 | 0.945 | 0.932 | 0.935 | 0.935 |

Group readings on epil's three cells that the default covers and the candidate misses, of 150. Group means are not scored
apart by the plan, so these count in G above.

| candidate | lost | true mean below 0.01 | truth below the interval | median true mean | median lower end, count scale |
|---|---|---|---|---|---|
| chi(4, 1.03) | 19 | 17 | 18 | 0.00019 | 0.0007 |
| chi(5, 0.9) | 25 | 20 | 22 | 0.00020 | 0.0013 |
| C4 at 0.6 | 42 | 35 | 40 | 0.00033 | 0.0032 |
| C4 at 0.75 | 49 | 36 | 47 | 0.00038 | 0.0089 |
| C5 at 1.3 | 35 | 31 | 33 | 0.00031 | 0.0015 |
| chi(3, 1.75) | 9 | 9 | 9 | 0.00009 | 0.0003 |
| chi(2.75, 2) | 9 | 9 | 9 | 0.00009 | 0.0003 |

On those cells, groups at or above 0.01 are covered 0.81 by the default and 0.82 by both priors. Groups below 0.01 are
covered 0.64 to 0.66 by the default, the figure depending on which replications are own runs, and 0.47 and 0.49 by the
two priors, against 0 to 0.35 for the five.

### The default-length check

This is step 5, what users would get. It ran for chi(2.75, 2) only: the check stops at the first prior that holds. Each
class's gap less the default's, at default lengths, must be below 0.02 at its 95 percent upper bound.

| class | gap less C0's (standard error) | 95 percent upper bound | status |
|---|---|---|---|
| F | +0.003 (0.001) | 0.005 | below the bar |
| G | +0.004 (0.004) | 0.011 | below the bar |
| R | +0.000 (0.005) | 0.009 | below the bar |
| P | -0.001 (0.001) | 0.001 | below the bar |

Coverage at default length, the default then chi(2.75, 2):

| class | narrow | middle | wide | small |
|---|---|---|---|---|
| F | 0.943, 0.941 | 0.918, 0.917 | 0.921, 0.924 | 0.894, 0.888 |
| G | 0.931, 0.931 | 0.840, 0.838 | 0.871, 0.875 | 0.809, 0.800 |
| R | 0.894, 0.897 | 0.868, 0.867 | 0.935, 0.950 | 0.941, 0.936 |
| P | 0.921, 0.920 | 0.888, 0.888 | 0.897, 0.898 | 0.897, 0.894 |
| F, rows below 0.01 | - | 0.617, 0.622 | 0.725, 0.711 | 0.726, 0.688 |

The check reads classes and not cells. Beside it, G on a few small-sample cells moves more: quakes-n100-3 from 0.850 to
0.775, epil-n30-3 from 0.483 to 0.400, roaches-n30-3 from 0.850 to 0.817. Each is 10 replications. Sampling takes about
the same time as under the default.

### What the two priors buy and give up

Bought: the prior check, which the default fails. On long runs, intervals on the log mean and on rate ratios score
lower, though less than the five's (0.82 and 0.83 under chi(2.75, 2) against 0.61 to 0.77); at default length much of
this goes, as the next section reads. As with the five most of this comes from one
design: on D4-3 F's score is 64.3 under the default, 33.9 under chi(2.75, 2) and 30.7 under chi(3, 1.75), against 6.7 to
24.0 for the five, its coverage 0.54, 0.60 and 0.62 against 0.63 to 0.77. Without D4-3 the ratios are 0.94 and 0.92, and
the median cell's is about 1. RMSE of the log mean falls on the small band from 1.77 to 1.27 and 1.22, and on the wide band
from 0.494 to 0.485 and 0.483, less than the five's.

Kept: narrow data, where the tried chi priors lost. D5, the null data, has interval width 0.49 under chi(3, 1.75) and
0.47 under chi(2.75, 2), against the default's 0.48.

Given up: G on epil's small and widened cells, by the amounts above, and R on the middle band by about 0.01. Near zero,
F's coverage below 0.01 falls by 0.03 to 0.04 on small data and by about 0.01 on wide on the long reading, and about the
same at default length. Against the five these are the smallest of those losses.

### What it cost and what a change would take

**Against the stop.** 15.2 in the stop's units in the harness's seconds, 15.4 in the workers', for a total of 109.7 to 109.9
of 110. In box core-hours, about 10.8. The reading scripts are not charged, as the first look's were not. No room is left.

**A fresh test.** Another block on chi(2.75, 2) alone would cost about 19 box core-hours on every band, 27 in the stop's
units; or about 5 and 7 without the wide band. Either passes the stop. Its boundary would be about 2.5, and it would be the
first reading of these priors on data not used to choose them. A block without the wide band would not retest epil-w3,
the largest gap.

**In the package, read and not built.** It is more than two numbers. Binary families and negative binomial share one
default line and one flag, so that line must learn the family, and binary families stay on chi(1.5, 2). That is about the
size of the plan's fixed-k change: some 30 lines of R, the help, the prior-defaults note and the tests. The engine takes the
prior's constants as data, so there is no engine change, and the shipped header carries no such prior, so no interface
change. `chi()` with no arguments stays chi(1.5, 2), which means that under negative binomial it would no longer equal the
default; the same prior on the sd is invchi(2.75, 1.5).

Every default negative binomial fit's draws move, and binary fits' do not. The exact gates are rerun, the negative
binomial mixing gate and its equivalence scenario are re-recorded on the reference build, and a design note is due. The
seeded-drift snapshots hold no negative binomial fit, so they do not move. The exact negative binomial target fixes k, and
the simulation-based calibration arm for this family does too, so no gate today samples the new prior's draw of k: such an
arm should be added.

### What was not checked

- **No fresh replications.** The two priors were picked after a rough look at these same outputs, so readings on the same
  replications can flatter them by an amount this reading cannot size; the next section's block was run to size it.
- **Own runs the rule requires that did not fit: 28 of 209,** the 8 on mite-bart, 19 on folds and one on mite-w3. These stay
  on the stored draws. Another 16 replications with a flat weight on chi(2.75, 2), whose weight against the default grows as
  a power of k, stay on the stored draws too.
- **chi(3, 1.75) has no default-length runs.** The check stopped at chi(2.75, 2).
- **Nothing was built.** Both builds used are the study's own, with the prior passed in.

## Added after a fresh block: chi(2.75, 2) on new replications

**The result.** The rule keeps the default. The prior chi(2.75, 2) was tested on 359 new replications, stated in advance in
two readings, and it had to pass both as well as the default-length check. It passes on all replications together, and it
passes the default-length check on all replications and on the new ones alone. On the new replications alone it passes
everything but two statistics. The interval score of group means is 1.020 against a limit of 1.05, with a standard error
of 0.022, so it sits under its limit by 1.37 standard errors where 1.96 were required. The coverage gap of rate ratios is
+0.0072 against a bar of 0.02, with a standard error of 0.0074, z 1.73.

So this is a failure to show the limits are met, not an estimate past them. Narrow data is unharmed, and no cell is worse
by the bar on either reading, though at 5 replications a cell that rule catches only losses that repeat in nearly every
replication. The interval score of group means is the statistic that was weakest at the first look, its noise comes from
the same few small-sample cells, and the group means on `epil`'s small and widened cells lose again on the new data. The
gap on rate ratios comes from one replication of one designed cell. The maintainer kept the default for now
([The outcome, and what a reopened study would need](#the-outcome-and-what-a-reopened-study-would-need)).

### What the search screened and found

The first look's candidates all failed. Every chi prior but one held the default's median k and raised the degrees of
freedom, which takes mass from near zero, as intended, and also from above k = 4, where k's posterior lies on null data.
Raising the scale as well removes mass near zero and adds it at the top. So the search looked along that direction before
any new fit: degrees of freedom 2 to 4 by 0.25 against scale 1.25 to 3.5 by 0.25, 90 settings. It read only the prior check
and draws stored from the first look, so it is a screen and not a reading of the rule. Eight settings pass the prior check:

| prior | spread limit, at most 1 | real-data floor, at least 1 | widest narrow-data ratio (limit 1.10) | cells worse by the bar | narrow F interval score |
|---|---|---|---|---|---|
| chi(2.25, 3) | 0.87 | 1.10 | 1.22 Traffic-bart | | 0.982 |
| chi(2.5, 2.25) | 0.99 | 1.33 | 1.07 Traffic-bart | Traffic-bart | 0.986 |
| chi(2.5, 2.5) | 0.77 | 1.01 | 1.11 Traffic-bart | Traffic-bart | 0.981 |
| chi(2.75, 2) | 0.82 | 1.20 | 1.04 Traffic-bart | | 0.994 |
| chi(3, 1.75) | 0.75 | 1.13 | 1.02 D5 | | 1.005 |
| chi(3.25, 1.5) | 0.79 | 1.19 | 1.07 D5 | epil-w3 | 1.022 |
| chi(3.5, 1.25) | 0.98 | 1.66 | 1.14 D5 | epil-n30-3, epil-w3 | 1.049 |
| chi(3.75, 1.25) | 0.69 | 1.18 | 1.14 D5 | epil-n30-3, epil-w3 | 1.047 |

All eight are tied on the rule's rank, within 0.003 of each other with rate ratios' designed truths corrected (0.002
before), against 0.005 for a tie. A larger scale loses on Traffic-bart, and a smaller one on null data's width and on
`epil-w3`. Two are not simple cases. chi(2.25, 3) escapes the flag only through noise: its weights keep as few as 4 effective
draws, so stored draws cannot read it, and its median k is twice the default's. chi(3, 1.75) ties with chi(2.75, 2) on
everything the rule reads, and the rule's last tie-break, the interval score on narrow data, gives chi(2.75, 2) the lead.
That pick was written into the plan before any new fit, along with the block and both readings below.

The prior's shape:

| | median k | P(k < 0.5) | P(k > 4) |
|---|---|---|---|
| chi(1.5, 2), the default | 1.91 | 0.080 | 0.085 |
| chi(2.75, 2) | 2.91 | 0.007 | 0.227 |

So it changes the bulk of the prior, not only its tails.

### What the block did

Five new replications on every cell that is not a fold, each cell continuing from where the first look stopped. On each, the
default's long run; chi(2.75, 2) read by reweighting that run on all its draws where the weights hold, and from its own run
where they do not and on every narrow replication; and default-length runs of both priors. The same bars and bounds as the
first look apply.

The two readings:

1. **All replications**: the first look's and the new ones together, at the spending boundary for the information reached.
   A statistic still pending does not pass. The information fraction is 0.686 and the boundary 2.47. No status depends on the
   small choices in computing it (2.46 to 2.47), nor on the corrected rate-ratio scores, which give 0.681 and 2.48; the
   smallest z is 3.03.
2. **The new replications alone**, free of how the prior was chosen: z against 1.96. The log score rests on the folds, which do
   not grow, so it is not part of this one.

The rule recommends the prior only if both readings and the default-length check pass.

### The readings

Candidate less default (standard error), z and status; the bar is 0.02 for gaps and 1.05 for interval score ratios. The last
column is the new replications less the first look on the same cells, how far the first look's reading moved.

| statistic | first look, boundary 3.46 | all replications, 2.47 | new alone, 1.96 | new less first look |
|---|---|---|---|---|
| F | +0.0002 (0.0008), 26.0, pass | +0.0007 (0.0006), 32.1, pass | +0.0012 (0.0010), 19.5, pass | +0.0010 (0.0012) |
| G | +0.0035 (0.0019), 8.6, pass | +0.0035 (0.0022), 7.6, pass | +0.0067 (0.0034), 3.9, pass | +0.0032 (0.0039) |
| R | -0.0063 (0.0027), 9.6, pass | -0.0058 (0.0046), 5.65, pass | **+0.0072 (0.0074), 1.73, not passed** | +0.0135 (0.0079) |
| P | -0.0005 (0.0010), 20.7, pass | +0.0005 (0.0008), 25.2, pass | -0.0001 (0.0007), 30.2, pass | +0.0005 (0.0012) |
| widest F | -0.0008 (0.0013), 16.5, pass | -0.0007 (0.0012), 16.9, pass | +0.0003 (0.0013), 14.6, pass | +0.0011 (0.0018) |
| widest G | +0.0081 (0.0045), 2.68, pending | +0.0039 (0.0026), 6.1, pass | -0.0000 (0.0020), 10.3, pass | -0.0081 (0.0049) |
| widest R | -0.0168 (0.0101), 3.63, pass | -0.0097 (0.0068), 4.3, pass | +0.0005 (0.0008), 23.1, pass | +0.0173 (0.0102) |
| widest P | +0.0001 (0.0009), 21.3, pass | -0.0012 (0.0008), 28.1, pass | -0.0019 (0.0008), 26.1, pass | -0.0020 (0.0013) |
| interval score F | 0.820 (0.019), 12.0, pass | 0.846 (0.017), 12.4, pass | 0.903 (0.026), 5.8, pass | +0.083 (0.032) |
| interval score G | 0.990 (0.022), 2.73, pending | 1.001 (0.016), 3.03, pass | **1.020 (0.022), 1.37, not passed** | +0.030 (0.031) |
| interval score R | 0.830 (0.015), 14.3, pass | 0.850 (0.014), 14.7, pass | 0.893 (0.023), 6.7, pass | +0.063 (0.028) |

The held-out log score, all replications: +0.0002 per row, upper bound 0.0007 against 0.01. Rank, the mean of F's and R's gaps:
-0.0030 at the first look, +0.0042 on the new replications.

The narrow-data guard, F's interval score ratio with its 95 percent interval, limit 1.10. It passes on every cell on each reading:

| cell | first look | all replications | new alone |
|---|---|---|---|
| D5 | 0.978 (0.974, 0.981) | 0.975 (0.971, 0.980) | 0.973 (0.963, 0.983) |
| D8 | 0.990 (0.985, 0.995) | 0.991 (0.989, 0.994) | 0.993 (0.991, 0.995) |
| Insurance | 0.982 (0.967, 0.997) | 0.987 (0.979, 0.994) | 0.992 (0.984, 1.000) |
| Insurance-bart | 0.981 (0.965, 0.996) | 0.986 (0.977, 0.994) | 0.990 (0.976, 1.005) |
| Traffic | 0.989 (0.974, 1.005) | 0.984 (0.977, 0.992) | 0.979 (0.973, 0.985) |
| Traffic-bart | 1.027 (1.004, 1.050) | 1.031 (1.019, 1.042) | 1.034 (1.016, 1.053) |
| quakes | 0.994 (0.989, 0.999) | 0.995 (0.992, 0.997) | 0.996 (0.992, 1.000) |
| quakes-bart | 0.994 (0.988, 1.000) | 0.995 (0.989, 1.001) | 0.996 (0.982, 1.011) |

**Cells worse by the bar: none on either reading**, and none among the shape-below-1 cells. A cell is worse when its lower bound
passes 0.02, by a t test on the cell's replications less one. At 5 replications a loss of the same size in each replication it
touches is flagged only if it touches at least 4 of the 5, so a clean result here is weak evidence; the largest differences stay
under the bound. The ones to look at:

| cell, statistic | first look | all replications | new alone |
|---|---|---|---|
| Traffic-bart, F | +0.026 (0.005), lower 0.012 | +0.022 (0.003), lower 0.015 | +0.019 (0.004), lower 0.008 |
| epil-w3, G | +0.105 (0.041), lower -0.009 | +0.050 (0.025), lower -0.007 | +0.000 (0.000), lower -0.000 |
| epil-n30-3, G | +0.052 (0.043), lower -0.046 | +0.057 (0.031), lower -0.010 | +0.133 (0.041), lower **0.0195** |
| epil-n100-3, G | +0.034 (0.022), lower -0.017 | +0.033 (0.018), lower -0.005 | +0.033 (0.033), lower -0.058 |
| roaches-w2, G | +0.000 (0.000), lower -0.000 | +0.000 (0.017), lower -0.038 | +0.095 (0.033), lower 0.003 |

`epil-n30-3`'s group means on the new replications miss the flag by 0.0005. `epil-w3`, the first look's largest gap and the
reason the wide band ran, does not repeat. Traffic-bart's F, close to the bar in the plan, stays under it.

### The two statistics that did not pass

Interval score of group means is a ratio of band-weighted averages of a log-scale score that charges 20 times a miss, and
small-sample cells set nearly all of its noise: half from `epil-n30-3`, 38 percent from D4-3, the design where the default's
intervals are worst, and 10 percent from `epil-n100-3`. The small band, a quarter of the weight by the plan, holds 68 percent of
the default's score on the new replications, though no real dataset in the study is both that small and that wide. By cell the median ratio
is 0.998, and only three of 67 cells stand above 1.05: `epil-n30-3` 1.17, `epil-n100-3` 1.14 and `epil-w3` 1.12. The same three
stood at 1.25, 1.27 and 1.12 at the first look. D4-3 pulls the other way, at 0.79. Read beside the rule, without that one
design, the ratio is 1.055 at the first look, 1.051 on all replications and 1.045 (0.019) on the new ones. So on all
replications the pass rests on D4-3, and the loss on `epil`'s group means is real and repeated: the rule's result is a
balance of two small sets of synthetic cells. Read beside the rule, on the new replications:

| reading | ratio (standard error) | z against 1.05 |
|---|---|---|
| as registered | 1.020 (0.022) | 1.37 |
| groups with true mean below 0.001 set apart | 0.983 (0.020) | 3.44 |
| below 0.01 set apart | 1.002 (0.015) | 3.18 |
| below 0.03 set apart | 1.018 (0.018) | 1.79 |
| below 0.1 set apart | 1.019 (0.018) | 1.73 |
| without `epil`'s cells | 0.954 (0.022) | 4.39 |
| without D4-3 | 1.045 (0.019) | 0.29 |
| without any designed cell | 1.050 (0.021) | 0.01 |
| without the seven cells whose spread passes 4 | 0.959 (0.020) | 4.51 |

That the loss sits on near-zero groups is only partly right: setting groups below 0.01 apart passes, but setting more apart
does not, and the readings do not move in order with the line, because D4-3's gain also sits on low groups. The last row
names a tension: the prior check calls a spread above 4 one where no count data have been seen, yet seven cells of the study
have a spread above 4, `epil`'s three behind this statistic among them. Either those truths are off-limits, and the rule's
keep dissolves, or they are not, and the check's reason for failing the default weakens.

To pass on the new replications the ratio had to be at most 1.008. With the standard error observed, a block of this size
passes this one statistic with a chance of about 0.64 if the two priors were equal and 0.28 if the true ratio were 1.02; the
bootstrap at 5 replications a cell understates its standard error by about a tenth, and with that corrected z is 1.23 and the
chances 0.54 and 0.24. The 95 percent interval, about 0.98 to 1.06, holds both 1.00 and 1.05.

The coverage gap of rate ratios, +0.0072 (0.0074), comes from one replication of D8, a narrow designed cell. There the
sample-average rate ratio's interval from chi(2.75, 2)'s own run misses its truth while the default's covers it, and R's
coverage on that replication falls from 0.93 to 0.42, half of R's score being that one interval. The same replication read by
reweighting the default's run covers it, at 0.93, so the miss is the own run's noise as much as the prior's. Before the
correction this gap read -0.0002 (0.0026) and passed.

### What the new replications say about the first look

The first look's readings were read after the prior was chosen on the same outputs, so the new replications are the way to look
for that. On the gaps, the bar the prior had to meet, the difference is within noise on every statistic, the largest being
widest R at +0.017 (0.010), R at +0.014 (0.008) and widest G at -0.008 (0.005).

The gain in interval score shrank: F from 0.820 to 0.903, a difference of 0.083 (0.032); R from 0.830 to 0.893, 0.063 (0.028);
G from 0.990 to 1.020, 0.030 (0.031). Most of that is one design, and on the default's side. D4-3's F score under the default
was 64 at the first look and 42 on the new replications, and its replications run from 10 to 97, so the change cannot be told
from noise; the default was not the prior chosen. Without D4-3 the differences are F +0.019 (0.014), G -0.011 (0.029), R
+0.013 (0.013), all within noise. So the new replications show no sign that choosing the prior on the first look's outputs
flattered it. They cannot rule it out either: they reuse the same cells and truths the search was read on, D4-3 being one
function from one seed, so they cure noise in the replications and not choice among the cells.

### Coverage by band

90 percent coverage on the long reading, the default then chi(2.75, 2), cells equal within a band.

| class | | narrow | middle | wide | small |
|---|---|---|---|---|---|
| F | first look | 0.947, 0.944 | 0.926, 0.923 | 0.931, 0.931 | 0.888, 0.891 |
| F | new | 0.937, 0.934 | 0.932, 0.931 | 0.930, 0.930 | 0.884, 0.886 |
| G | first look | 0.931, 0.931 | 0.848, 0.843 | 0.891, 0.887 | 0.838, 0.832 |
| G | new | 0.925, 0.946 | 0.833, 0.819 | 0.877, 0.877 | 0.791, 0.790 |
| R | first look | 0.897, 0.896 | 0.892, 0.881 | 0.961, 0.949 | 0.971, 0.959 |
| R | new | 0.904, 0.877 | 0.941, 0.949 | 0.972, 0.973 | 0.989, 0.978 |
| P | first look | 0.917, 0.917 | 0.900, 0.900 | 0.895, 0.895 | 0.900, 0.891 |
| P | new | 0.925, 0.925 | 0.899, 0.899 | 0.898, 0.899 | 0.868, 0.867 |
| F, rows below 0.01 | first look | - | 0.654, 0.653 | 0.848, 0.839 | 0.886, 0.854 |
| F, rows below 0.01 | new | - | 0.717, 0.721 | 0.877, 0.867 | 0.876, 0.856 |

The first look's P here leaves out the folds, which the new replications do not have. R's new narrow reading is D8's one
replication above. The cost near zero repeats, at 0.02 to 0.03 on
small data and 0.01 on wide; rows below 0.01 are not read by the rule.

### The default-length check

What users would get, each class's gap less the default's, which must be below 0.02 at its 95 percent upper bound.

| class | first look | all replications | new alone |
|---|---|---|---|
| F | +0.0029 (0.0012), upper 0.0052 | +0.0021 (0.0010), upper 0.0041 | +0.0005 (0.0018), upper 0.0040 |
| G | +0.0040 (0.0037), upper 0.0113 | +0.0016 (0.0035), upper 0.0084 | -0.0025 (0.0053), upper 0.0079 |
| R | +0.0001 (0.0046), upper 0.0091 | +0.0013 (0.0039), upper 0.0091 | -0.0010 (0.0034), upper 0.0057 |
| P | -0.0010 (0.0009), upper 0.0009 | -0.0002 (0.0008), upper 0.0014 | -0.0018 (0.0016), upper 0.0014 |

It passes on both, on the 222 of the 359 new replications that ran. The check reads classes and not cells. Beside it, on all
replications, R on D7 is +0.040 (0.005), past the bar at the cell's bound, while D7's three new replications alone read
-0.015 (0.011). The three small-sample cells that moved at the first look do not repeat: `epil-n30-3`'s G goes from 0.278 to
0.333 on the new replications, `quakes-n100-3`'s stays at 0.583 and `roaches-n30-3`'s at 0.833, each on three. Sampling takes
the same time under either prior, within 1 percent.

### What the prior buys at default length

The rule reads interval scores on the long runs only; the default-length check reads coverage gaps. Read beside the rule, the
interval scores at default length, chi(2.75, 2) over the default (standard error):

| statistic | first look | new | all replications | long runs, all replications |
|---|---|---|---|---|
| interval score F | 0.965 (0.011) | 0.980 (0.013) | 0.968 (0.009) | 0.846 (0.017) |
| interval score G | 1.042 (0.019) | 1.005 (0.026) | 1.031 (0.015) | 1.001 (0.016) |
| interval score R | 0.943 (0.015) | 0.955 (0.009) | 0.945 (0.012) | 0.850 (0.014) |

Most of the long-run gain does not reach what users run: k has not mixed at default length, and D4-3's F score under the default
is 15 and 13 there against 64 and 42 on long runs. Read by the rule's own statistics at default length, chi(2.75, 2) would not pass
on the first look, the new replications or all together: on all replications G's interval score is 1.031 (z 1.20), widest R's gap has z 0.97, and D7's rate ratios are
worse by the bar.

### At the spreads of real data

Read beside the rule, by kind of cell: the 22 real datasets at their own spread (as fitted and with a BART-fitted truth), the 19
widened ones, the 20 small samples and the 7 designed cells. Median cell interval score ratio, chi(2.75, 2) over the default, then
coverage less the default's, all replications:

| class | own spread, long | own spread, default length | widened, long | small, long | designed, long |
|---|---|---|---|---|---|
| F | 0.995, -0.003 | 0.996, -0.002 | 0.995, +0.000 | 0.983, +0.003 | 0.996, -0.001 |
| G | 0.996, -0.000 | 1.000, -0.004 | 0.998, -0.003 | 0.999, -0.004 | 0.997, -0.005 |
| R | 0.994, -0.005 | 0.990, -0.004 | 0.984, -0.008 | 0.961, -0.012 | 0.994, -0.002 |

At their own spread the two priors are practically the same. On those cells' new replications the median shift in a row's
posterior mean of the log mean is 0.004 between the priors, 0.02 of its posterior sd, against 0.014 between a default-length and a
long run under the default; for group means, 0.002 against 0.004. Run length moves estimates more than this prior does.

### What the rule does

Reading 1 passes. Reading 2 does not, on the interval score of group means and the coverage gap of rate ratios. The
default-length check passes. The rule recommends chi(2.75, 2) only if all three pass, so **it keeps the default**, chi(1.5, 2), which still fails the prior
check. This was the study's last look: the rule funds no more runs.

### What could not be run

- **`mite-bart`, all five new replications,** cut in the plan before the block so that the rest would fit, `epil-w3` being
  costlier but kept.
- **`mite-w3`'s fifth new replication:** the stop left no room for it.
- **137 of the 359 default-length replications.** They were held in reserve, but a few long runs on widened cells went well past
  what the first look had priced, one `mite-w3` replication by more than three times. What remained ran cell by cell, first new
  replication everywhere and so on, so every cell is read: 60 cells at 3 of 5, the 9 cheapest small-sample cells at 4,
  and `quakes`, `quakes-bart` and `quakes-w3` at 2. The rest would have cost about 0.9 more units.
- **The folds,** which never grow. The log score is the first look's.
- **Any prior but chi(2.75, 2) on new data,** chi(3, 1.75) included, as the plan states.

The weights held throughout, keeping at least 0.66 of the draws' effective share, with k-hat at most 0.24. This prior's weight
against the default grows without bound as a power of k, which matters where k's posterior reaches high; no replication outside
the narrow band had more than 5 percent of its draws above k = 4, and every narrow replication has its own run. 99 own runs ran, 40
on narrow data for that reason and 59 of 319 elsewhere because the rule asks. On the narrow band the reweighted and own readings
agree: F's interval score within 0.5 percent and its coverage within 0.01 on every cell.

### Fixed k

The study tried only two fixed values, k = 1.5 and 2. They fail the prior check, at 0.19 and about 0 on the real-data floor,
against a floor of 1, so no fit read them. When the maintainer asked whether fixed k had been checked, lower values were read the same way, on
default-length runs of 5 replications a cell. k = 1, 1.1 and 1.2 pass the prior check, at 0.72, 0.24 and 0.04 percent on the
spread limit and 4.25, 2.27 and 1.32 on the real-data floor. They fail the narrow-data guard by far. On null data, F's interval
score is 1.89 to 1.96 times the default's against a limit of 1.10, and intervals are about 1.8 times as wide. Seven of the eight
narrow cells fail at each value. A fixed k small enough for the prior check gives up what the default does on data with little to
find, so fixed k is not a way out.

### Cost

27.10 units, 19.08 box core-hours, 68,700 worker-seconds against a stop of 68,690, so 10 seconds past it. In the harness's seconds
it is 26.91. That brings the study to 137.0 of the stop of 137.

## Point estimates and effects, read after the fact

The rule reads coverage and interval scores. A user of `family = "nbinom"` reads a point estimate and its interval for a
quantity of the question at hand: a mean count, a group's mean, a treatment's rate ratio. Everything in this section was
chosen after the results were known and decides nothing.

The quantities: the log mean at a row and group means, by RMSE and bias; for rate ratios, the per-row effect (its RMSE,
sometimes called PEHE), effects averaged within groups, and the sample-average rate ratio, by its error, its bias and its
interval's width. Each is chi(2.75, 2) over the default: the ratio of the cells' means within a kind, then the median cell's
ratio. Bias is the cell's mean error over replications, taken in absolute value. The long runs are the new replications'; the
default-length runs are every replication's.

| quantity | own spread | widened | small | designed |
|---|---|---|---|---|
| long runs | | | | |
| log mean, RMSE | 0.997, 0.996 | 0.982, 0.996 | 0.826, 0.977 | 0.994, 0.997 |
| log mean, bias | 0.959, 0.972 | 0.952, 0.959 | 0.723, 0.871 | 0.949, 0.969 |
| group means, RMSE | 0.999, 0.998 | 1.018, 0.998 | 0.968, 0.995 | 1.002, 1.003 |
| group means, bias | 0.972, 0.985 | 1.036, 0.977 | 0.898, 0.996 | 1.003, 1.002 |
| per-row effect, RMSE | 0.994, 1.000 | 0.995, 1.001 | 0.935, 0.976 | 0.999, 0.997 |
| effects within groups, RMSE | 0.995, 0.996 | 0.999, 1.005 | 0.928, 0.980 | 0.999, 0.997 |
| sample-average effect, error | 0.998, 1.001 | 1.001, 0.999 | 0.972, 0.993 | 1.000, 1.004 |
| sample-average effect, bias | 1.008, 1.005 | 1.002, 0.990 | 1.024, 0.994 | 1.010, 1.012 |
| sample-average effect, width | 0.983, 0.985 | 0.970, 0.986 | 0.836, 0.935 | 0.993, 0.995 |
| default length | | | | |
| log mean, RMSE | 0.997, 0.999 | 0.989, 0.991 | 0.941, 0.990 | 0.999, 0.999 |
| log mean, bias | 0.991, 1.032 | 0.948, 0.931 | 0.773, 0.852 | 0.934, 0.999 |
| group means, RMSE | 0.992, 0.999 | 1.004, 0.993 | 1.021, 1.005 | 1.002, 1.005 |
| group means, bias | 0.917, 0.981 | 1.013, 1.025 | 1.102, 1.023 | 0.955, 0.960 |
| per-row effect, RMSE | 0.994, 0.999 | 0.989, 0.995 | 0.942, 0.988 | 0.990, 0.992 |
| effects within groups, RMSE | 1.016, 1.002 | 1.009, 0.989 | 0.930, 0.971 | 0.972, 0.986 |
| sample-average effect, error | 0.993, 0.994 | 1.047, 0.959 | 0.937, 0.977 | 1.021, 0.995 |
| sample-average effect, bias | 0.942, 0.990 | 1.053, 0.994 | 1.005, 1.017 | 0.955, 1.006 |
| sample-average effect, width | 0.981, 0.990 | 0.995, 0.985 | 0.921, 0.954 | 0.981, 0.999 |

Off the small samples every ratio of RMSE, error and width is within 0.96 to 1.05, and bias moves by under 10 percent
either way from values small to begin with: the default's mean absolute bias of the log mean at its own spread is 0.016,
against an RMSE of 0.25. On the small samples,
which hold D4-2 and D4-3, chi(2.75, 2) gains on the log mean and on effects, by about 17 percent on the log mean's RMSE on long
runs and 6 percent at default length; the median small cell gains 2 percent and 1 percent.

On the real outcomes, held-out rows read on the first look's long runs: a row's posterior mean of the log mean moves by a
median of 0.0006 to 0.013 between the priors, by dataset, against 0.004 to 0.26 between a default-length and a long run
under the default (`mite` 0.26, `crimtab` 0.18). A group's mean moves by 0.0001 to 0.0055, at most 0.012. The new count's
coverage moves by at most 0.002, and the held-out log score over a whole dataset by -1.25 to +0.43 nats.

Two results hold under both priors alike, and so say nothing about the choice between them. 13 of the 40 cells with rate
ratios have a true effect that varies across rows: the designed cells and five truths fitted by BART. On D7 the true per-row
effect has sd 0.18 and the per-row effect's RMSE is 0.22 under each prior, so a single effect for all rows would have done
better. On D3 the sample-average log rate ratio is 0.16 low under each prior, a rate ratio about 15 percent low, with an
interval 0.38 wide.

The quantities themselves were not chosen for each dataset's question. Group means are taken over the levels of the first
column with 2 to 10 values and over quartiles of the first numeric column, and the rate ratio contrasts the first binary
column. On `roaches`, `epil` and `Traffic`, and arguably `solder` and `Arabidopsis`, these are what an analyst would ask; on
`quakes`, `crimtab`, `mite`, `Insurance` and `aids` they are arbitrary, and five of the 11 real outcomes have no rate ratio.

## The outcome, and what a reopened study would need

Under the rule the default stays: k drawn under chi(1.5, 2). It fails the prior check, and no candidate that passes it
passed the fits, the sixth only by not showing two limits were met on the new replications. The maintainer kept the default
for now and left the question open. Nothing is built.

What the study does settle:

- **Fixed k is not a way out.** Values small enough for the prior check widen intervals on data with little to find by
  about 1.8 times.
- **Floors on k cost group means near zero,** and tighter chi priors cost width on null data; neither is a candidate.
- **At the spreads of real data the choice barely matters** to intervals, point estimates or effects, and run length
  matters more than this prior.

What it does not settle is whether chi(2.75, 2), or another prior that passes the prior check, should replace the default.
The evidence for and against it sits on synthetic cells far from the real datasets, on a reading users do not run, at a
number of replications that cannot decide the statistic that kept the default. A study that reopens the question would need:

1. **Datasets people meet,** chosen by a stated rule and not for being wide, each read at its own spread with a truth fitted
   as it is and one fitted by BART. Here the small and widened cells carry most of the noise and most of the weight, and no real
   dataset in the study is both small and wide.
2. **Estimands natural to each dataset's question,** stated per dataset: group means on the count scale, the sample-average
   rate ratio and count difference, per-row and group effects; scored by RMSE and bias as well as by coverage and interval
   score.
3. **A reading of what users run, decided first.** Most of what chi(2.75, 2) buys on long runs is gone at default length, and
   the rule read interval scores on long runs only.
4. **Enough replications to decide.** At 5 a cell, the cell rule flags only losses that repeat in at least 4 of the 5, and the
   class statistic that kept the default had a chance of only 0.54 to 0.64 of passing if the two priors were equal.
5. **A criterion argued as others argue for a default prior.** The prior check's limit, at most 1 in 100 generated datasets
   with a spread above 4, has no published anchor, and the limits on interval scores, 1.05 and 1.10, and the trade between
   classes were set without a derivation. Comparable practice scales a forest's prior from the data, as Murray (2021) does from
   a near-maximum count, or keeps heavy tails, and judges a default by predictive accuracy and estimand error. A scale set
   from the data was tried here only through a GLM's spread; a near-maximum-count setting was never scored.
6. **Real outcomes read for more than the new count.** They enter only the coverage of a new count and the held-out log
   score, and the cells draw independent counts at one shape, so clustering, zero inflation beyond the negative binomial,
   shape varying with the predictors and serial dependence are untested.

At the fresh block's measured prices, 12 more datasets with two truths each at 5 replications would take about 4 box
core-hours, with default-length runs and own runs where the weights fail.

## Disclosures

- **Unsettled cells stayed in the band averages.** The plan reports a cell with more than 20 percent of its
  replications unsettled apart from its band's average. The harness that made the decision kept `mite-bart` and
  mite's folds in. Recomputed without them, on the same scores and boundary, all five verdicts and the outcome are
  the same.
- **C6 went mostly unscored because of a harness defect.** A k-hat the tail fit cannot return was counted as a failed
  check; see [The shape's prior](#the-shapes-prior). The five fitted candidates' k-hat is never missing, so their own-run
  choices and the decision are not touched.
- **Rate ratios' truth on seven designed cells was stored as zero.** On D1, D2, D3, D4-2, D4-3, D6 and D8 the stored
  truth for R read the binary covariate's column in place of the effect's size, so it was zero while the counts carried the
  effect. Every rate-ratio reading here is scored against the corrected truth, from the stored draws with no new fit, and no
  verdict moves. The information fractions and boundaries computed from the corrected scores move slightly, to 0.381 and
  3.45 at the first look and 0.681 and 2.48 on all replications; the tables keep the registered 3.46 and 2.47, and no status
  differs at either. On the first look C4 at 0.75 is no longer flagged on D4-3's rate ratios, though still on its group means;
  on the fresh block the coverage gap of rate ratios on the new replications alone, which passed before, is no longer shown
  below its bar; under [Sensitivity](#sensitivity) chi(7, 0.75) now fails.
- These are the three departures found between the plan and the code. G scoring its near-zero groups with the rest is
  not one: the plan sets such rows apart in F and P only. Its effect is in [Sensitivity](#sensitivity).
- **The price.** Look 1 cost about twice its sizing, for the reasons under [Cost](#cost). The stop was not reached,
  and this look's result does not depend on it: every candidate fails whatever the stop reads.
  At measured prices a block of 5 replications on every band costs 63.4 in the table's units against the 37
  pre-registered, of which 15.5 was left; no block was called for then. The fresh block ran later, on a stop the
  maintainer raised.
- **Pending parts.** G and the widest cells' G are pending for several candidates. They do not matter here, since
  each of those candidates has already failed on another part.
- **The rescoring for this note.** The tables here come from read-only rescoring of the same outputs. It reproduces
  the decision's statistics exactly with the stored truths, and writes no decision.

## Readings by cell

The default's coverage and gap on the long reading, its coverage at default length, then each candidate's gap less
the default's on the long reading; `*` marks a
cell worse by the bar. `m` is the cell's replications; a dash means no replication the candidate could be
paired on.

F, the log mean at a row:

| cell | band | m | C0 cover | C0 gap | C0 default cover | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 |
|---|---|---|---|---|---|---|---|---|---|---|
| D5 | narrow | 5 | 0.988 | 0.089 | 0.986 | -0.005 | -0.005 | +0.000 | +0.000 | +0.012 |
| D8 | narrow | 5 | 0.952 | 0.053 | 0.946 | +0.005 | +0.004 | +0.000 | +0.000 | -0.002 |
| Insurance | narrow | 5 | 1.000 | 0.101 | 1.000 | -0.000 | +0.000 | +0.000 | +0.000 | -0.019 |
| Insurance-bart | narrow | 5 | 0.988 | 0.088 | 0.988 | +0.000 | +0.003 | +0.000 | +0.000 | -0.022 |
| quakes | narrow | 5 | 0.960 | 0.062 | 0.959 | -0.002 | -0.003 | +0.000 | +0.000 | -0.001 |
| quakes-bart | narrow | 5 | 0.943 | 0.044 | 0.937 | +0.002 | +0.004 | +0.000 | +0.000 | -0.004 |
| Traffic | narrow | 5 | 0.966 | 0.067 | 0.968 | +0.001 | -0.001 | +0.000 | +0.000 | -0.067 |
| Traffic-bart | narrow | 5 | 0.778 | 0.121 | 0.762 | -0.044 | -0.059 | +0.000 | +0.000 | +0.183 * |
| Arabidopsis | middle | 5 | 0.927 | 0.028 | 0.929 | +0.001 | +0.002 | +0.000 | +0.000 | +0.002 |
| Arabidopsis-bart | middle | 5 | 0.847 | 0.052 | 0.820 | -0.004 | -0.004 | +0.000 | +0.000 | -0.003 |
| D1 | middle | 5 | 0.913 | 0.013 | 0.899 | +0.001 | -0.000 | +0.000 | +0.000 | -0.001 |
| D2 | middle | 5 | 0.946 | 0.048 | 0.937 | -0.002 | -0.004 | +0.000 | -0.000 | +0.000 |
| D3 | middle | 5 | 0.916 | 0.019 | 0.910 | +0.001 | +0.001 | +0.000 | +0.000 | +0.000 |
| D6 | middle | 5 | 0.947 | 0.053 | 0.928 | -0.002 | +0.002 | -0.002 | -0.004 | +0.002 |
| D7 | middle | 5 | 0.936 | 0.038 | 0.934 | +0.000 | -0.000 | +0.000 | +0.000 | +0.000 |
| epil | middle | 5 | 0.961 | 0.062 | 0.961 | +0.003 | +0.003 | +0.000 | +0.000 | +0.001 |
| epil-bart | middle | 5 | 0.870 | 0.026 | 0.869 | -0.005 | -0.005 | +0.000 | +0.000 | -0.004 |
| grouseticks | middle | 5 | 0.932 | 0.029 | 0.925 | -0.000 | -0.000 | +0.000 | +0.000 | -0.000 |
| grouseticks-bart | middle | 5 | 0.908 | 0.000 | 0.889 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches | middle | 5 | 0.907 | 0.000 | 0.911 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches-bart | middle | 5 | 0.976 | 0.077 | 0.960 | -0.001 | -0.001 | +0.000 | +0.000 | +0.000 |
| solder | middle | 5 | 0.980 | 0.081 | 0.979 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids | wide | 5 | 0.930 | 0.031 | 0.920 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids-bart | wide | 5 | 0.910 | 0.000 | 0.898 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids-w2 | wide | 5 | 0.935 | 0.036 | 0.931 | -0.001 | -0.002 | +0.000 | +0.000 | -0.002 |
| aids-w3 | wide | 5 | 0.914 | 0.005 | 0.885 | -0.005 | -0.005 | +0.000 | -0.005 | -0.005 |
| Arabidopsis-w2 | wide | 5 | 0.959 | 0.059 | 0.951 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Arabidopsis-w3 | wide | 5 | 0.969 | 0.070 | 0.966 | -0.000 | -0.004 | +0.000 | +0.003 | -0.000 |
| crimtab | wide | 10 | 0.958 | 0.059 | 0.950 | -0.002 | -0.002 | +0.000 | -0.002 | -0.001 |
| crimtab-bart | wide | 10 | 0.941 | 0.044 | 0.932 | -0.001 | -0.001 | -0.003 | -0.012 | -0.002 |
| crimtab-w3 | wide | 10 | 0.942 | 0.044 | 0.928 | -0.001 | -0.001 | -0.005 | -0.015 | -0.003 |
| epil-w2 | wide | 5 | 0.932 | 0.033 | 0.919 | -0.003 | -0.005 | +0.000 | +0.000 | +0.000 |
| epil-w3 | wide | 5 | 0.884 | 0.000 | 0.925 | +0.000 | +0.012 | +0.000 | +0.035 | +0.022 |
| grouseticks-w2 | wide | 5 | 0.959 | 0.060 | 0.955 | -0.002 | -0.002 | +0.000 | +0.000 | +0.000 |
| grouseticks-w3 | wide | 5 | 0.967 | 0.067 | 0.950 | -0.002 | -0.009 | +0.000 | +0.000 | +0.000 |
| mite | wide | 5 | 0.947 | 0.050 | 0.942 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| mite-bart | wide | 5 | 0.895 | 0.000 | 0.852 | +0.000 | +0.000 | +0.003 | +0.001 | +0.000 |
| mite-w2 | wide | 5 | 0.965 | 0.067 | 0.955 | +0.001 | -0.003 | +0.000 | +0.000 | -0.000 |
| mite-w3 | wide | 5 | 0.849 | 0.047 | 0.864 | -0.007 | -0.008 | -0.047 | -0.047 | -0.047 |
| quakes-w2 | wide | 5 | 0.951 | 0.052 | 0.929 | +0.000 | -0.000 | +0.000 | +0.000 | +0.000 |
| quakes-w3 | wide | 5 | 0.942 | 0.045 | 0.921 | -0.001 | +0.005 | +0.000 | +0.000 | +0.000 |
| roaches-w2 | wide | 5 | 0.905 | 0.000 | 0.891 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches-w3 | wide | 5 | 0.943 | 0.044 | 0.944 | -0.002 | -0.004 | -0.002 | -0.007 | -0.000 |
| solder-bart | wide | 5 | 0.866 | 0.032 | 0.853 | +0.002 | +0.002 | +0.000 | -0.000 | -0.000 |
| solder-w2 | wide | 5 | 0.958 | 0.059 | 0.953 | -0.001 | -0.001 | +0.000 | +0.000 | -0.000 |
| solder-w3 | wide | 5 | 0.933 | 0.035 | 0.924 | +0.005 | +0.005 | -0.001 | +0.003 | +0.002 |
| Traffic-w2 | wide | 5 | 0.923 | 0.023 | 0.916 | -0.002 | -0.003 | +0.000 | -0.001 | -0.001 |
| Traffic-w3 | wide | 5 | 0.926 | 0.027 | 0.896 | +0.002 | -0.005 | -0.002 | -0.003 | -0.001 |
| aids-n100-2 | small | 10 | 0.898 | 0.000 | 0.891 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids-n100-3 | small | 10 | 0.908 | 0.000 | 0.895 | +0.009 | +0.015 | +0.010 | +0.007 | +0.012 |
| aids-n30-2 | small | 10 | 0.880 | 0.006 | 0.870 | +0.008 | +0.008 | -0.000 | +0.008 | -0.006 |
| aids-n30-3 | small | 10 | 0.880 | 0.000 | 0.876 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| D4-2 | small | 10 | 0.915 | 0.014 | 0.929 | +0.022 | +0.026 | +0.018 | +0.012 | +0.020 |
| D4-3 | small | 10 | 0.540 | 0.356 | 0.660 | -0.099 | -0.108 | -0.117 | -0.227 | -0.088 |
| epil-n100-2 | small | 10 | 0.925 | 0.025 | 0.924 | -0.007 | +0.000 | +0.001 | +0.002 | +0.001 |
| epil-n100-3 | small | 10 | 0.876 | 0.000 | 0.930 | +0.000 | +0.009 | +0.034 | +0.028 | +0.040 |
| epil-n30-2 | small | 10 | 0.893 | 0.000 | 0.883 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| epil-n30-3 | small | 10 | 0.899 | 0.000 | 0.899 | +0.000 | +0.000 | +0.000 | +0.030 | +0.000 |
| Insurance-w2 | small | 10 | 0.997 | 0.097 | 0.995 | +0.000 | +0.000 | +0.000 | +0.000 | -0.002 |
| Insurance-w3 | small | 10 | 0.989 | 0.089 | 0.986 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| quakes-n100-2 | small | 10 | 0.911 | 0.010 | 0.907 | +0.002 | +0.002 | +0.000 | +0.000 | -0.000 |
| quakes-n100-3 | small | 10 | 0.912 | 0.010 | 0.914 | +0.005 | +0.002 | +0.002 | +0.009 | +0.004 |
| quakes-n30-2 | small | 10 | 0.943 | 0.044 | 0.940 | -0.010 | -0.013 | +0.000 | -0.007 | +0.000 |
| quakes-n30-3 | small | 10 | 0.866 | 0.029 | 0.856 | +0.009 | -0.000 | -0.001 | -0.021 | +0.000 |
| roaches-n100-2 | small | 10 | 0.890 | 0.000 | 0.881 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches-n100-3 | small | 10 | 0.892 | 0.000 | 0.888 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches-n30-2 | small | 10 | 0.863 | 0.027 | 0.887 | +0.014 | +0.022 | -0.016 | +0.022 | +0.008 |
| roaches-n30-3 | small | 10 | 0.876 | 0.018 | 0.866 | -0.003 | -0.018 | +0.003 | +0.021 | -0.002 |

G, the mean count over a group:

| cell | band | m | C0 cover | C0 gap | C0 default cover | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 |
|---|---|---|---|---|---|---|---|---|---|---|
| D5 | narrow | 5 | 0.967 | 0.067 | 0.967 | -0.000 | -0.055 | +0.000 | +0.000 | +0.034 |
| D8 | narrow | 5 | 0.867 | 0.023 | 0.867 | +0.036 | +0.036 | +0.000 | +0.000 | -0.000 |
| Insurance | narrow | 5 | 0.950 | 0.037 | 0.950 | -0.000 | +0.063 | +0.000 | +0.000 | +0.000 |
| Insurance-bart | narrow | 5 | 1.000 | 0.101 | 1.000 | -0.000 | +0.000 | +0.000 | +0.000 | -0.063 |
| quakes | narrow | 5 | 0.850 | 0.009 | 0.850 | -0.000 | -0.001 | +0.000 | +0.000 | -0.000 |
| quakes-bart | narrow | 5 | 0.950 | 0.038 | 0.950 | +0.000 | +0.000 | +0.000 | +0.000 | -0.000 |
| Traffic | narrow | 5 | 0.933 | 0.011 | 0.933 | +0.000 | +0.000 | +0.000 | +0.000 | -0.011 |
| Traffic-bart | narrow | 5 | 0.933 | 0.011 | 0.933 | +0.001 | +0.000 | +0.000 | +0.000 | -0.011 |
| Arabidopsis | middle | 5 | 0.933 | 0.000 | 0.867 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Arabidopsis-bart | middle | 5 | 0.800 | 0.087 | 0.800 | -0.000 | -0.000 | +0.000 | +0.000 | -0.000 |
| D1 | middle | 5 | 0.833 | 0.034 | 0.833 | +0.000 | -0.000 | +0.000 | -0.000 | +0.040 |
| D2 | middle | 5 | 0.900 | 0.000 | 0.867 | +0.022 | +0.022 | +0.000 | +0.000 | +0.000 |
| D3 | middle | 5 | 0.733 | 0.145 | 0.667 | -0.000 | -0.000 | +0.000 | +0.000 | -0.000 |
| D6 | middle | 5 | 0.700 | 0.192 | 0.733 | -0.000 | -0.097 | -0.000 | -0.066 | -0.000 |
| D7 | middle | 5 | 0.433 | 0.466 | 0.467 | -0.000 | -0.000 | +0.000 | +0.000 | -0.000 |
| epil | middle | 5 | 0.933 | 0.000 | 0.967 | +0.067 | +0.067 | +0.000 | +0.000 | +0.067 |
| epil-bart | middle | 5 | 0.867 | 0.000 | 0.867 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| grouseticks | middle | 5 | 0.886 | 0.000 | 0.886 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| grouseticks-bart | middle | 5 | 0.914 | 0.000 | 0.914 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches | middle | 5 | 0.967 | 0.067 | 0.933 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches-bart | middle | 5 | 0.967 | 0.067 | 0.967 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| solder | middle | 5 | 1.000 | 0.100 | 1.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids | wide | 5 | 0.875 | 0.000 | 0.875 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids-bart | wide | 5 | 0.950 | 0.037 | 0.950 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids-w2 | wide | 5 | 0.875 | 0.000 | 0.875 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids-w3 | wide | 5 | 0.875 | 0.000 | 0.900 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Arabidopsis-w2 | wide | 5 | 0.933 | 0.000 | 0.933 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Arabidopsis-w3 | wide | 5 | 1.000 | 0.101 | 1.000 | -0.000 | +0.000 | +0.000 | +0.001 | -0.000 |
| crimtab | wide | 10 | 1.000 | 0.101 | 1.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| crimtab-bart | wide | 10 | 0.900 | 0.000 | 0.900 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| crimtab-w3 | wide | 10 | 0.925 | 0.000 | 0.875 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| epil-w2 | wide | 5 | 0.767 | 0.104 | 0.767 | +0.000 | -0.030 | -0.000 | -0.000 | -0.000 |
| epil-w3 | wide | 5 | 0.767 | 0.128 | 0.567 | +0.070 | +0.105 | +0.202 | +0.235 * | +0.238 * |
| grouseticks-w2 | wide | 5 | 0.914 | 0.000 | 0.829 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| grouseticks-w3 | wide | 5 | 0.914 | 0.000 | 0.914 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| mite | wide | 5 | 0.891 | 0.000 | 0.891 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| mite-bart | wide | 5 | 0.909 | 0.000 | 0.927 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| mite-w2 | wide | 5 | 0.891 | 0.000 | 0.855 | +0.000 | +0.013 | +0.000 | +0.000 | +0.000 |
| mite-w3 | wide | 5 | 0.909 | 0.000 | 0.800 | +0.013 | +0.000 | +0.078 | +0.134 | +0.005 |
| quakes-w2 | wide | 5 | 0.800 | 0.078 | 0.750 | -0.000 | -0.000 | +0.000 | +0.000 | -0.000 |
| quakes-w3 | wide | 5 | 0.900 | 0.000 | 0.900 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches-w2 | wide | 5 | 0.967 | 0.067 | 0.933 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches-w3 | wide | 5 | 1.000 | 0.101 | 0.900 | -0.000 | +0.000 | +0.000 | -0.032 | +0.000 |
| solder-bart | wide | 5 | 0.867 | 0.000 | 0.933 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| solder-w2 | wide | 5 | 0.800 | 0.088 | 0.800 | +0.000 | +0.000 | +0.000 | -0.000 | -0.000 |
| solder-w3 | wide | 5 | 0.933 | 0.000 | 0.933 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Traffic-w2 | wide | 5 | 0.700 | 0.199 | 0.733 | -0.000 | -0.000 | -0.000 | -0.000 | -0.000 |
| Traffic-w3 | wide | 5 | 0.900 | 0.000 | 0.900 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids-n100-2 | small | 10 | 0.863 | 0.000 | 0.875 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids-n100-3 | small | 10 | 0.900 | 0.000 | 0.875 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids-n30-2 | small | 10 | 0.762 | 0.136 | 0.750 | -0.027 | -0.027 | -0.013 | -0.027 | +0.013 |
| aids-n30-3 | small | 10 | 0.900 | 0.000 | 0.863 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| D4-2 | small | 10 | 0.800 | 0.095 | 0.750 | +0.034 | +0.016 | -0.029 | +0.017 | +0.016 |
| D4-3 | small | 10 | 0.833 | 0.030 | 0.817 | -0.030 | -0.030 | +0.079 | +0.352 * | +0.062 |
| epil-n100-2 | small | 10 | 0.750 | 0.147 | 0.733 | +0.018 | +0.018 | -0.000 | -0.000 | -0.000 |
| epil-n100-3 | small | 10 | 0.783 | 0.115 | 0.633 | +0.034 | +0.101 * | +0.199 * | +0.250 * | +0.118 |
| epil-n30-2 | small | 10 | 0.767 | 0.132 | 0.733 | +0.017 | -0.000 | +0.017 | -0.018 | +0.017 |
| epil-n30-3 | small | 10 | 0.733 | 0.163 | 0.483 | +0.219 * | +0.252 * | +0.318 * | +0.385 * | +0.303 * |
| Insurance-w2 | small | 10 | 1.000 | 0.100 | 1.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Insurance-w3 | small | 10 | 1.000 | 0.100 | 1.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| quakes-n100-2 | small | 10 | 0.775 | 0.124 | 0.775 | -0.000 | -0.000 | -0.000 | +0.000 | -0.000 |
| quakes-n100-3 | small | 10 | 0.800 | 0.087 | 0.850 | -0.048 | -0.048 | -0.049 | +0.026 | -0.000 |
| quakes-n30-2 | small | 10 | 0.875 | 0.000 | 0.850 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| quakes-n30-3 | small | 10 | 0.850 | 0.042 | 0.850 | +0.030 | +0.019 | +0.031 | +0.102 | +0.031 |
| roaches-n100-2 | small | 10 | 0.850 | 0.044 | 0.833 | -0.000 | +0.019 | -0.000 | -0.000 | -0.000 |
| roaches-n100-3 | small | 10 | 0.817 | 0.080 | 0.800 | -0.000 | -0.000 | -0.001 | -0.001 | -0.000 |
| roaches-n30-2 | small | 10 | 0.867 | 0.000 | 0.850 | +0.053 | +0.052 | +0.000 | +0.053 | +0.000 |
| roaches-n30-3 | small | 10 | 0.833 | 0.053 | 0.850 | -0.032 | -0.053 | -0.034 | -0.022 | -0.000 |

R, rate ratios of a binary covariate:

| cell | band | m | C0 cover | C0 gap | C0 default cover | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 |
|---|---|---|---|---|---|---|---|---|---|---|
| D5 | narrow | 5 | 0.755 | 0.092 | 0.749 | -0.018 | -0.010 | +0.000 | +0.000 | +0.009 |
| D8 | narrow | 5 | 0.988 | 0.089 | 0.983 | +0.003 | +0.006 | +0.000 | +0.000 | -0.005 |
| Traffic | narrow | 5 | 0.992 | 0.092 | 0.994 | -0.004 | -0.006 | +0.000 | +0.000 | -0.092 |
| Traffic-bart | narrow | 5 | 0.853 | 0.000 | 0.849 | +0.000 | +0.000 | +0.000 | +0.000 | +0.162 |
| Arabidopsis | middle | 5 | 1.000 | 0.101 | 0.995 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Arabidopsis-bart | middle | 5 | 0.802 | 0.055 | 0.789 | -0.014 | -0.014 | +0.000 | +0.000 | -0.012 |
| D1 | middle | 5 | 0.966 | 0.068 | 0.954 | +0.000 | +0.001 | +0.000 | +0.000 | +0.001 |
| D2 | middle | 5 | 0.720 | 0.155 | 0.710 | +0.007 | +0.007 | -0.000 | +0.001 | -0.000 |
| D3 | middle | 5 | 0.770 | 0.088 | 0.655 | -0.000 | +0.001 | +0.000 | +0.000 | -0.000 |
| D6 | middle | 5 | 0.968 | 0.071 | 0.960 | -0.001 | -0.000 | -0.001 | +0.000 | +0.007 |
| D7 | middle | 5 | 0.744 | 0.114 | 0.742 | +0.000 | +0.000 | +0.000 | +0.000 | -0.001 |
| epil | middle | 5 | 1.000 | 0.101 | 1.000 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| epil-bart | middle | 5 | 0.796 | 0.076 | 0.789 | -0.003 | -0.003 | +0.000 | +0.000 | -0.009 |
| roaches | middle | 5 | 0.950 | 0.041 | 0.945 | +0.000 | +0.000 | +0.000 | +0.000 | +0.001 |
| roaches-bart | middle | 5 | 0.998 | 0.099 | 0.888 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| solder | middle | 5 | 0.990 | 0.090 | 0.990 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Arabidopsis-w2 | wide | 5 | 0.787 | 0.055 | 0.878 | -0.000 | -0.000 | +0.000 | -0.000 | -0.000 |
| Arabidopsis-w3 | wide | 5 | 0.991 | 0.092 | 0.993 | +0.004 | +0.004 | +0.002 | -0.001 | +0.002 |
| epil-w2 | wide | 5 | 1.000 | 0.102 | 0.899 | -0.000 | -0.102 | +0.000 | +0.000 | +0.000 |
| epil-w3 | wide | 5 | 1.000 | 0.101 | 0.790 | -0.101 | -0.101 | -0.101 | -0.084 | -0.101 |
| mite | wide | 5 | 0.998 | 0.100 | 0.989 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| mite-bart | wide | 5 | 0.978 | 0.082 | 0.954 | -0.009 | -0.008 | -0.004 | -0.009 | -0.005 |
| mite-w2 | wide | 5 | 1.000 | 0.101 | 1.000 | +0.000 | +0.001 | +0.000 | +0.001 | +0.000 |
| mite-w3 | wide | 5 | 0.999 | 0.101 | 0.892 | -0.000 | +0.001 | +0.001 | +0.000 | +0.000 |
| roaches-w2 | wide | 5 | 1.000 | 0.101 | 0.997 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches-w3 | wide | 5 | 1.000 | 0.102 | 1.000 | -0.000 | +0.000 | +0.000 | +0.001 | +0.000 |
| solder-bart | wide | 5 | 0.918 | 0.020 | 0.907 | +0.000 | +0.000 | +0.000 | +0.000 | +0.001 |
| solder-w2 | wide | 5 | 0.983 | 0.083 | 0.980 | -0.001 | +0.000 | +0.000 | +0.000 | +0.000 |
| solder-w3 | wide | 5 | 0.986 | 0.087 | 0.984 | +0.001 | -0.000 | +0.001 | +0.001 | -0.000 |
| Traffic-w2 | wide | 5 | 0.779 | 0.064 | 0.765 | +0.001 | -0.005 | -0.000 | -0.001 | -0.001 |
| Traffic-w3 | wide | 5 | 1.000 | 0.102 | 1.000 | +0.000 | -0.000 | +0.000 | -0.000 | +0.000 |
| D4-2 | small | 10 | 0.990 | 0.091 | 0.982 | -0.091 | -0.091 | -0.091 | -0.024 | -0.091 |
| D4-3 | small | 10 | 1.000 | 0.103 | 0.955 | -0.001 | -0.010 | -0.103 | -0.059 | -0.066 |
| epil-n100-2 | small | 10 | 1.000 | 0.102 | 0.999 | +0.000 | -0.000 | +0.000 | +0.000 | +0.000 |
| epil-n100-3 | small | 10 | 1.000 | 0.102 | 0.949 | -0.064 | -0.065 | -0.069 | -0.102 | -0.064 |
| epil-n30-2 | small | 10 | 0.987 | 0.087 | 0.983 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| epil-n30-3 | small | 10 | 1.000 | 0.101 | 0.883 | +0.000 | -0.068 | -0.101 | -0.041 | -0.101 |
| roaches-n100-2 | small | 10 | 0.906 | 0.000 | 0.907 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches-n100-3 | small | 10 | 0.933 | 0.000 | 0.884 | +0.000 | +0.000 | +0.000 | +0.086 | +0.000 |
| roaches-n30-2 | small | 10 | 0.913 | 0.000 | 0.903 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches-n30-3 | small | 10 | 0.983 | 0.084 | 0.963 | +0.000 | +0.000 | +0.000 | -0.006 | +0.000 |

P, a new count at a held-out row, on truths and on the folds. D4-3 has no row here: in all 10 replications all 12 of
its held-out rows have a true mean below 0.01 (the largest 0.0015), so they are scored apart.

| cell | band | m | C0 cover | C0 gap | C0 default cover | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 |
|---|---|---|---|---|---|---|---|---|---|---|
| D5 | narrow | 5 | 0.884 | 0.013 | 0.880 | +0.001 | -0.003 | +0.000 | +0.000 | +0.013 |
| D8 | narrow | 5 | 0.887 | 0.000 | 0.888 | +0.009 | +0.007 | +0.000 | +0.000 | +0.000 |
| Insurance | narrow | 5 | 0.954 | 0.054 | 0.954 | +0.015 | +0.015 | +0.000 | +0.000 | +0.000 |
| Insurance-bart | narrow | 5 | 0.969 | 0.069 | 0.969 | +0.000 | +0.000 | +0.000 | +0.000 | -0.015 |
| Insurance-folds | narrow | 5 | 0.985 | 0.085 | 0.985 | +0.000 | -0.015 | +0.000 | +0.000 | -0.017 |
| quakes | narrow | 5 | 0.892 | 0.000 | 0.896 | +0.006 | +0.009 | +0.000 | +0.000 | +0.000 |
| quakes-bart | narrow | 5 | 0.920 | 0.019 | 0.919 | +0.000 | +0.000 | +0.000 | +0.000 | -0.002 |
| quakes-folds | narrow | 5 | 0.910 | 0.010 | 0.908 | -0.003 | -0.002 | +0.000 | +0.000 | -0.003 |
| Traffic | narrow | 5 | 0.914 | 0.000 | 0.914 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Traffic-bart | narrow | 5 | 0.914 | 0.000 | 0.914 | +0.000 | +0.000 | +0.000 | +0.000 | +0.010 |
| Traffic-folds | narrow | 5 | 0.902 | 0.000 | 0.902 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Arabidopsis | middle | 5 | 0.904 | 0.000 | 0.906 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Arabidopsis-bart | middle | 5 | 0.893 | 0.000 | 0.890 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Arabidopsis-folds | middle | 5 | 0.818 | 0.082 | 0.821 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| D1 | middle | 5 | 0.914 | 0.010 | 0.920 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| D2 | middle | 5 | 0.888 | 0.011 | 0.890 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| D3 | middle | 5 | 0.873 | 0.026 | 0.877 | +0.004 | +0.004 | +0.000 | +0.000 | +0.000 |
| D6 | middle | 5 | 0.897 | 0.000 | 0.897 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| D7 | middle | 5 | 0.891 | 0.006 | 0.891 | +0.000 | +0.000 | +0.000 | +0.000 | -0.003 |
| epil | middle | 5 | 0.906 | 0.000 | 0.902 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| epil-bart | middle | 5 | 0.915 | 0.006 | 0.915 | -0.006 | -0.006 | +0.000 | +0.000 | -0.006 |
| epil-folds | middle | 5 | 0.881 | 0.018 | 0.877 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| grouseticks | middle | 5 | 0.904 | 0.000 | 0.906 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| grouseticks-bart | middle | 5 | 0.916 | 0.016 | 0.919 | +0.000 | +0.000 | +0.000 | +0.000 | -0.002 |
| grouseticks-folds | middle | 5 | 0.908 | 0.006 | 0.916 | -0.002 | +0.001 | -0.006 | +0.002 | -0.006 |
| roaches | middle | 5 | 0.919 | 0.000 | 0.919 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches-bart | middle | 5 | 0.873 | 0.024 | 0.877 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches-folds | middle | 5 | 0.771 | 0.129 | 0.767 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| solder | middle | 5 | 0.908 | 0.000 | 0.906 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| solder-folds | middle | 5 | 0.893 | 0.004 | 0.887 | +0.001 | +0.003 | +0.000 | +0.000 | +0.000 |
| aids | wide | 5 | 0.908 | 0.000 | 0.910 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids-bart | wide | 5 | 0.914 | 0.010 | 0.910 | -0.004 | -0.004 | +0.000 | +0.000 | +0.000 |
| aids-folds | wide | 5 | 0.892 | 0.000 | 0.890 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids-w2 | wide | 5 | 0.897 | 0.000 | 0.895 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids-w3 | wide | 5 | 0.883 | 0.000 | 0.888 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Arabidopsis-w2 | wide | 5 | 0.901 | 0.000 | 0.896 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Arabidopsis-w3 | wide | 5 | 0.902 | 0.000 | 0.902 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| crimtab | wide | 10 | 0.907 | 0.006 | 0.909 | +0.001 | +0.002 | +0.002 | +0.000 | +0.002 |
| crimtab-bart | wide | 10 | 0.901 | 0.000 | 0.904 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| crimtab-folds | wide | 5 | 0.925 | 0.025 | 0.925 | +0.001 | -0.000 | +0.003 | +0.002 | -0.002 |
| crimtab-w3 | wide | 10 | 0.898 | 0.000 | 0.898 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| epil-w2 | wide | 5 | 0.945 | 0.045 | 0.945 | +0.000 | -0.009 | +0.000 | +0.000 | +0.000 |
| epil-w3 | wide | 5 | 0.789 | 0.111 | 0.789 | +0.000 | +0.000 | -0.025 | -0.040 | +0.000 |
| grouseticks-w2 | wide | 5 | 0.909 | 0.001 | 0.909 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| grouseticks-w3 | wide | 5 | 0.896 | 0.000 | 0.899 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| mite | wide | 5 | 0.901 | 0.000 | 0.902 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| mite-bart | wide | 5 | 0.894 | 0.005 | 0.896 | +0.001 | +0.001 | -0.005 | -0.001 | -0.001 |
| mite-folds | wide | 5 | 0.918 | 0.017 | 0.914 | +0.002 | +0.002 | +0.002 | +0.001 | +0.001 |
| mite-w2 | wide | 5 | 0.899 | 0.000 | 0.899 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| mite-w3 | wide | 5 | 0.894 | 0.000 | 0.911 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| quakes-w2 | wide | 5 | 0.912 | 0.010 | 0.908 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| quakes-w3 | wide | 5 | 0.892 | 0.008 | 0.892 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches-w2 | wide | 5 | 0.877 | 0.023 | 0.877 | -0.009 | -0.004 | +0.000 | -0.004 | +0.000 |
| roaches-w3 | wide | 5 | 0.896 | 0.000 | 0.904 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| solder-bart | wide | 5 | 0.884 | 0.015 | 0.884 | +0.000 | +0.000 | +0.000 | +0.000 | -0.001 |
| solder-w2 | wide | 5 | 0.903 | 0.000 | 0.902 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| solder-w3 | wide | 5 | 0.892 | 0.000 | 0.884 | +0.000 | +0.000 | +0.000 | +0.000 | +0.007 |
| Traffic-w2 | wide | 5 | 0.870 | 0.030 | 0.870 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Traffic-w3 | wide | 5 | 0.903 | 0.000 | 0.903 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids-n100-2 | small | 10 | 0.840 | 0.057 | 0.840 | +0.006 | +0.006 | +0.000 | +0.000 | -0.007 |
| aids-n100-3 | small | 10 | 0.911 | 0.000 | 0.895 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids-n30-2 | small | 10 | 0.883 | 0.000 | 0.883 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| aids-n30-3 | small | 10 | 0.963 | 0.063 | 0.947 | -0.063 | -0.063 | -0.018 | -0.063 | -0.018 |
| D4-2 | small | 10 | 0.925 | 0.022 | 0.925 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| epil-n100-2 | small | 10 | 0.925 | 0.018 | 0.915 | +0.000 | +0.000 | +0.000 | -0.018 | +0.000 |
| epil-n100-3 | small | 10 | 0.885 | 0.000 | 0.899 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| epil-n30-2 | small | 10 | 0.867 | 0.004 | 0.900 | +0.000 | -0.004 | +0.000 | +0.000 | +0.000 |
| epil-n30-3 | small | 10 | 0.917 | 0.000 | 0.917 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Insurance-w2 | small | 10 | 0.977 | 0.077 | 0.977 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Insurance-w3 | small | 10 | 0.977 | 0.077 | 0.977 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| quakes-n100-2 | small | 10 | 0.835 | 0.065 | 0.830 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| quakes-n100-3 | small | 10 | 0.916 | 0.000 | 0.926 | +0.010 | +0.010 | +0.010 | +0.000 | +0.010 |
| quakes-n30-2 | small | 10 | 0.967 | 0.067 | 0.950 | -0.067 | -0.067 | +0.000 | -0.017 | +0.000 |
| quakes-n30-3 | small | 10 | 0.922 | 0.000 | 0.897 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches-n100-2 | small | 10 | 0.860 | 0.039 | 0.865 | -0.006 | -0.006 | +0.000 | +0.000 | +0.000 |
| roaches-n100-3 | small | 10 | 0.845 | 0.055 | 0.865 | -0.005 | -0.005 | -0.005 | -0.010 | -0.005 |
| roaches-n30-2 | small | 10 | 0.833 | 0.064 | 0.800 | +0.034 | +0.034 | +0.016 | +0.034 | +0.034 |
| roaches-n30-3 | small | 10 | 0.850 | 0.033 | 0.833 | +0.028 | +0.049 | +0.049 | +0.033 | +0.028 |

## Sensitivity

This section is not an option offered after the fact. It shows what the criterion's other readings, shown to the
maintainer with stage 1's table, would have recommended on these results.

| reading of Part A | who passes A | what the rule gives |
|---|---|---|
| as written, A3 at the widest fit, 1 in 100 | none | the default kept |
| the widest fit, 1 in 200 | chi(4, 1.03), C4 at 0.6, C5 at 1.3 | all three fail step 2 above: the default kept |
| no room condition, A1 and A2 alone | eleven priors, fixed k included | on what ran, none passes; see below |
| the check reported only, a change needing a gain of 0.02 | all | no candidate gains 0.02 on F or R: the default kept |

Under the last two readings other priors would have had their own runs and did not. Read from weights alone,
where weights hold, which is where the prior matters least and so flatters them:

| prior | replications on truths, of 485 | cells | what step 2 reads |
|---|---|---|---|
| chi(3, 1.25) | 393 | 71 | pending: G's interval score 1.09 |
| chi(7, 0.75) | 175 | 42 | fails: R's interval score 1.10 |
| chi(10, 0.62) | 100 | 33 | fails: R's interval score 1.08 |
| chi(1.5, 3) | 439 | 65 | passes, gaps 0 and ratios 1.00, no narrow cell read |
| C4, k at least 1 | 193 | 43 | no reading on any of the widest cells |
| half-normal, median 1.57 | 353 | 61 | pending: G's interval score 1.18 |
| C7 | 251 | 53 | pending: G's interval score 1.05 |

Fixed k was not run. chi(1.5, 3) fails Part A's A1 at 4.62 percent under every reading but the last, and under the
last it gains nothing.

### What decides C4 at 0.6: found after scoring

This was found after the scores were seen. It is not offered as an option, and it cannot re-decide
anything by itself.

C4 at 0.6 fails only on groups whose true mean count is below about 0.001. The plan scores rows below 0.01 apart in F
and P, because how far each prior lets the log mean fall decides their coverage. It has no such rule for G, and the
code follows the plan. This table runs the decision's own code with G's groups below a line of true mean count set
apart, the boundary recomputed at each line from the corrected rate-ratio scores:

| line | information fraction | boundary | outcome | C4 at 0.6 |
|---|---|---|---|---|
| none, as registered | 0.381 | 3.45 | the default | fails |
| 0.001 | 0.381 | 3.45 | unresolved, the default or C4 at 0.6 | pending, widest R at z 3.44 |
| 0.003 | 0.382 | 3.44 | unresolved, the default or C4 at 0.6 | pending, widest R at z 3.44 |
| 0.01 | 0.389 | 3.41 | C4 at 0.6 | passes |
| 0.03 | 0.403 | 3.34 | C4 at 0.6 | passes |
| 0.1 | 0.403 | 3.34 | C4 at 0.6 | passes |

At the 0.01 line the pass is narrow: widest R's estimate is -0.016, in the candidate's favor, at z 3.44 against 3.41. G's
interval score ratio goes from 1.25 to 0.91. Both chi priors fail at every line on the width guard, and C5 on F in
Traffic-bart. C4 at 0.75 fails at the two lowest lines, on group means in D4-3 and `epil-n30-3`, and is pending from 0.01
up, on widest R at z 1.92; the outcome there is C4 at 0.6 under either resolution.

Near-zero means are still paid for. Under C4 at 0.6, F's coverage on rows below 0.01 falls from 0.89 to 0.68 on small
data and from 0.85 to 0.77 on wide, and G's on epil's groups below 0.01 from 0.66 to 0.02. Those rows are scored apart,
not ignored.

Reading G as F and P are read would be a deviation from the pre-registration, which only the maintainer can rule. Step
5 would still be unread: C4 at 0.6's default-length runs have not been run.

Dropping epil's widened and small-sample cells altogether, a coarser cut, leaves C4 at 0.6 pending and the look
unresolved between the default and C4 at 0.6, and the other four still fail. Those cells were fixed before any fit.

## Appendix: evidence

- **The pre-registration and its amendments,** in the study's plan,
  [nbinom-default-study.md](../plans/nbinom-default-study.md): the criterion and
  [The decision rule](../plans/nbinom-default-study.md#the-decision-rule),
  [Box readiness](../plans/nbinom-default-study.md#box-readiness),
  [Amendment before the fresh block](../plans/nbinom-default-study.md#amendment-before-the-fresh-block), which fixed
  chi(2.75, 2), the block and both readings before any fit on them, and
  [The defect in R's truth, and its rescoring](../plans/nbinom-default-study.md#the-defect-in-rs-truth-and-its-rescoring).
- **The maintainer's rulings,** in the decision ledger: dec-B412 and dec-B431 on the earlier prior predictive study,
  dec-B432 on this study and its rule, dec-B502 on reading A3 at the middle of a dataset's fits, dec-B556 on the search
  and the fresh block, and dec-B558 on keeping the default with the question open. dec-B559 sets the standard for a
  future study that changes a default prior.
- **The code,** the study's harness,
  [benchmarks/R/nbinom-default-study.R](../../benchmarks/R/nbinom-default-study.R), which generates the cells, runs the
  fits and scores them; every table here is read with it.
- **The stored outputs:** each run's 500 stored draws and readings, kept on the run machine. Every reading here, the
  rescored rate ratios and the reads after the fact among them, is read from those outputs, with no fit beyond the
  study's runs.
- **Reference:** J. S. Murray (2021). Log-linear Bayesian additive regression trees for multinomial logistic and count
  regression models. *Journal of the American Statistical Association* 116(534), 756-769.
