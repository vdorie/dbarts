# The negative binomial default prior: a study with its criterion fixed first

Status: FINDING, 2026-10-31 (dec-B432, dec-B502). The redesigned study of the default prior of
`family = "nbinom"`, run as pre-registered in [nbinom-default-study.md](../plans/nbinom-default-study.md). Its
first look decided: the rule keeps the default, k drawn under chi(1.5, 2). Awaiting the opus review against the
pre-registration, then the maintainer's ruling; nothing is built until then.

## Summary

The rule recommends keeping the default. The default fails Part A: 4.0 to 8.4 percent of its prior datasets
span a spread of the log mean above 4, against 1 percent allowed. Five candidates pass Part A, and all five fail
Part B at the first look. So what passing A would cost in fits is measured, and it is more than the rule
allows: each candidate loses coverage or interval score somewhere the rule forbids.

- **chi(4, 1.03)** and **chi(5, 0.9)** widen intervals on narrow and null data past the 10 percent guard (null
  data 21 and 26 percent), and lose coverage of group means on epil's small samples.
- **C4 at 0.6** and **C4 at 0.75**, the default with k kept above a floor, leave narrow data exactly as they
  are, and make the interval score of group means 25 and 39 percent worse.
- **C5 at 1.3**, the half-normal, under-covers one narrow cell (0.60 against the default's 0.78) and loses
  coverage of group means on epil's widened cells.

The losses center on cells whose truth is wider than the prior now allows. On epil widened to a spread of 3, a
spread over all rows of 4.43, past Part A's line of 4, the default's 90 percent intervals for group means
already cover only 0.73 to 0.78. The candidates' narrower intervals cover 0.35 to 0.75. On average they gain
elsewhere: F's and R's interval scores fall 22 to 39 percent, most of it on a few wide small-sample cells. The
rule does not trade a loss on one class or cell against a gain on another.

The first look used 94.5 of the 110 core-hours allowed in the sizing table's units. No guard run and no
confirmation block is called for.

## The question and the rule

The question (dec-B432): is there a better default prior on the forest's spread for negative binomial fits,
judged by a criterion fixed before any fit is scored? The criterion has two parts. Part A reads the prior alone.
Part B reads fits: coverage of 90 percent intervals on truths with known values, the interval score, and
held-out log score. A candidate that passes A replaces a default that fails A only if it loses nothing on B. The
rule, its bars and its bounds are in [The decision rule](../plans/nbinom-default-study.md#the-decision-rule) and
[Box readiness](../plans/nbinom-default-study.md#box-readiness), as pre-registered.

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
crimtab (dec-B502).

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
candidate to check, so no guard runs. The reading is final at this look.

Each candidate's gap less the default's (standard error) and status. The gap is the distance from the expected
coverage on either side, corrected for noise. The bar is 0.02.

| statistic | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 |
|---|---|---|---|---|---|
| F | -0.002 (0.001) pass | -0.003 (0.001) pass | -0.001 (0.001) pass | -0.002 (0.001) pass | +0.002 (0.002) pass |
| G | +0.007 (0.002) pass | +0.006 (0.003) pass | +0.009 (0.002) pass | +0.016 (0.002) pending | +0.010 (0.003) pending |
| R | -0.003 (0.003) pass | -0.009 (0.003) pass | -0.008 (0.002) pass | +0.002 (0.004) pass | -0.001 (0.015) pending |
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
| interval score R | 0.78 (0.018) pass | 0.74 (0.019) pass | 0.67 (0.026) pass | 0.65 (0.027) pass | 0.70 (0.025) pass |

The log score is C0's less the candidate's per row, the bar 0.01 at its 95 percent upper bound; all pass. Rank is
the mean of F's and R's gaps less C0's.

| row | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 |
|---|---|---|---|---|---|
| log score, C0 less candidate (95 upper) | +0.0030 (+0.0041) | +0.0041 (+0.0053) | -0.0003 (+0.0005) | +0.0004 (+0.0012) | +0.0008 (+0.0039) |
| rank | -0.0026 | -0.0058 | -0.0049 | -0.0001 | +0.0002 |

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
| C4 at 0.75 | R | D4-3 | +0.203 | 0.068 | 1.000 | 0.593 |
| C5 at 1.3 | G | epil-n30-3 | +0.303 | 0.065 | 0.733 | 0.433 |
| C5 at 1.3 | G | epil-w3 | +0.238 | 0.041 | 0.767 | 0.533 |
| C5 at 1.3 | F | Traffic-bart | +0.183 | 0.032 | 0.778 | 0.595 |

No shape-below-1 cell is worse by the bar for any candidate.

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
68 cells in the reviewer's constructed case, and 0.015 none. That does not touch these flags: the smallest flagged
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
| narrow | 0.898 | 0.906 | 0.904 | 0.898 | 0.898 | 0.872 | 0.899 |
| middle | 0.859 | 0.850 | 0.859 | 0.858 | 0.859 | 0.859 | 0.898 |
| wide | 0.961 | 0.954 | 0.948 | 0.954 | 0.947 | 0.954 | 0.898 |
| small | 0.972 | 0.967 | 0.962 | 0.945 | 0.906 | 0.955 | 0.899 |

P:

| band | C0 | chi(4, 1.03) | chi(5, 0.9) | C4 at 0.6 | C4 at 0.75 | C5 at 1.3 | expected |
|---|---|---|---|---|---|---|---|
| narrow | 0.921 | 0.921 | 0.920 | 0.921 | 0.921 | 0.916 | 0.900 |
| middle | 0.888 | 0.888 | 0.888 | 0.888 | 0.888 | 0.888 | 0.900 |
| wide | 0.897 | 0.897 | 0.896 | 0.898 | 0.899 | 0.897 | 0.900 |
| small | 0.900 | 0.892 | 0.892 | 0.897 | 0.895 | 0.897 | 0.900 |

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
| R | narrow | 0.898 | 0.069 | 0.895 | 0.071 | 0.35 | 0.34 |
| R | middle | 0.859 | 0.098 | 0.837 | 0.103 | 0.75 | 0.72 |
| R | wide | 0.961 | 0.086 | 0.935 | 0.062 | 1.85 | 1.50 |
| R | small | 0.972 | 0.068 | 0.942 | 0.041 | 6.58 | 3.59 |
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
small bands, and on the rows below 0.01 it falls from 0.85 to 0.73 on wide data and from 0.89 to 0.73 on small.

## What keeping the default costs

The default keeps a prior whose top tail spends its mass where no count data have been seen: 4.0 to 8.4 percent
of prior datasets past a spread of 4, by design, against the 1 percent of Part A. A change of default was expected
to buy a coherent prior predictive, not better intervals; that part stays unbought.

It also keeps very wide intervals on a few wide small-sample cells. The interval score of F on D4-3 (n = 60,
spread 3) is 64.3 under the default against 9.0 under C4 at 0.6, and 10.5 against 4.5 on epil-n100-3. These few
cells produce most of the fall in F's and R's interval score ratios: over the small band F's score is 7.18 under
the default and 3.22 to 4.57 under the candidates, and over the wide band 1.90 against 1.64 to 1.85. RMSE of the
log mean moves with it (wide band 0.494 against 0.398 to 0.476).

What the default does not keep is good coverage of group means where the truth is widest: 0.73 to 0.78 on epil's
widened cells, and 0.84 over the small band against an expected 0.90. Every candidate is worse there, not better.

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

**Mixing, for the collapsed move's plan.** Sweeps per effective draw of k on the scratch build, pooled over 4
chains, at each run's stop: median 28 on narrow, 89 on middle, 239 on wide, 157 on small, and 146 on the folds
with a 90th percentile of 1,519. Probes of today's sampler measured 440 to 620 on narrow data and 1,140 to 8,000
where the spread is 1 or more. The rule stops on effective draws and not on R-hat: k's split R-hat at the stop has
median 1.035, 90th percentile 1.09 and largest 1.55, and the mean count's largest is 1.09. On as-fitted real data
the median is 18 to 159, mite the largest. Trees on mite's real outcomes and BART-fitted truth mix slowly enough
that the rows, not k, set the length.

## The shape's prior

C6 is a separate decision and is read where its per-draw weights hold: 133 of 540 replications, 121 of them on
truths, in 31 cells, one narrow replication and 7 small among them. There its gap differences are -0.001 to 0, and
its interval score ratios 0.98 to 1.00. Elsewhere it is unscored. Its verdict is unfinished for that reason and says
nothing in favor of a change; it would take its own runs to read.

## Cost

Look 1 in box core-hours, against the pre-registered sizing at the calibrated speed. In the sizing table's units,
in which the stop is set, it comes to 94.5 of 110.

| component | sized | measured |
|---|---|---|
| the default's long runs | 12.67 | 21.18 |
| own runs | 16.41 | 42.97 |
| readings from all draws | 1.83 | 1.44 |
| the default's default readings | 1.55 | 0.96 |
| all | 32.46 | 66.55 |

| group | sized | measured |
|---|---|---|
| narrow | 3.31 | 3.73 |
| middle | 1.90 | 3.39 |
| wide | 18.45 | 32.87 |
| small | 5.28 | 7.88 |
| shape below 1 | 0.92 | 1.97 |
| folds | 2.54 | 16.92 |

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

## Disclosures

- **Unsettled cells stayed in the band averages.** The plan reports a cell with more than 20 percent of its
  replications unsettled apart from its band's average. The harness that made the decision kept `mite-bart` and
  mite's folds in. Recomputed without them, on the same scores and boundary, all five verdicts and the outcome are
  the same. It is the one departure this reading found between the plan and the decision's code.
- **The price.** Look 1 cost about twice its sizing, for the reasons under [Cost](#cost). The stop was not reached,
  and this look's result does not depend on it: every candidate fails whatever the stop reads.
  At measured prices a block of 5 replications on every band costs 63.4 in the table's units against the 37
  pre-registered, of which 15.5 was left; no block is called for.
- **Pending parts.** G and the widest cells' G are pending for several candidates. They do not matter here, since
  each of those candidates has already failed on another part.
- **The rescoring for this note.** The tables here come from a read-only rescoring of the same outputs. It reproduces
  the decision's statistics exactly and writes no decision or queue.

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
| D8 | narrow | 5 | 0.991 | 0.092 | 0.990 | -0.005 | -0.008 | +0.000 | +0.000 | -0.002 |
| Traffic | narrow | 5 | 0.992 | 0.092 | 0.994 | -0.004 | -0.006 | +0.000 | +0.000 | -0.092 |
| Traffic-bart | narrow | 5 | 0.853 | 0.000 | 0.849 | +0.000 | +0.000 | +0.000 | +0.000 | +0.162 |
| Arabidopsis | middle | 5 | 1.000 | 0.101 | 0.995 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| Arabidopsis-bart | middle | 5 | 0.802 | 0.055 | 0.789 | -0.014 | -0.014 | +0.000 | +0.000 | -0.012 |
| D1 | middle | 5 | 0.717 | 0.173 | 0.712 | +0.001 | +0.000 | +0.000 | -0.000 | -0.000 |
| D2 | middle | 5 | 0.739 | 0.115 | 0.734 | +0.002 | +0.003 | -0.000 | +0.000 | -0.000 |
| D3 | middle | 5 | 0.885 | 0.000 | 0.881 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| D6 | middle | 5 | 0.681 | 0.211 | 0.583 | +0.114 | -0.007 | +0.005 | -0.007 | +0.002 |
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
| D4-2 | small | 10 | 1.000 | 0.101 | 1.000 | +0.000 | -0.000 | -0.000 | -0.000 | +0.000 |
| D4-3 | small | 10 | 1.000 | 0.103 | 0.949 | -0.001 | -0.001 | -0.103 | +0.203 * | -0.001 |
| epil-n100-2 | small | 10 | 1.000 | 0.102 | 0.999 | +0.000 | -0.000 | +0.000 | +0.000 | +0.000 |
| epil-n100-3 | small | 10 | 1.000 | 0.102 | 0.949 | -0.064 | -0.065 | -0.069 | -0.102 | -0.064 |
| epil-n30-2 | small | 10 | 0.987 | 0.087 | 0.983 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| epil-n30-3 | small | 10 | 1.000 | 0.101 | 0.883 | +0.000 | -0.068 | -0.101 | -0.041 | -0.101 |
| roaches-n100-2 | small | 10 | 0.906 | 0.000 | 0.907 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches-n100-3 | small | 10 | 0.933 | 0.000 | 0.884 | +0.000 | +0.000 | +0.000 | +0.086 | +0.000 |
| roaches-n30-2 | small | 10 | 0.913 | 0.000 | 0.903 | +0.000 | +0.000 | +0.000 | +0.000 | +0.000 |
| roaches-n30-3 | small | 10 | 0.983 | 0.084 | 0.963 | +0.000 | +0.000 | +0.000 | -0.006 | +0.000 |

P, a new count at a held-out row, on truths and on the folds:

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
| chi(7, 0.75) | 175 | 42 | pending: R |
| chi(10, 0.62) | 100 | 33 | fails: R's interval score 1.07 |
| chi(1.5, 3) | 439 | 65 | passes, gaps 0 and ratios 1.00, no narrow cell read |
| C4, k at least 1 | 193 | 43 | no reading on any of the widest cells |
| half-normal, median 1.57 | 353 | 61 | pending: G's interval score 1.18 |
| C7 | 251 | 53 | pending: G's interval score 1.05 |

Fixed k was not run. chi(1.5, 3) fails Part A's A1 at 4.62 percent under every reading but the last, and under the
last it gains nothing.

How much rests on one family: without epil's widened and small-sample cells, C4 at 0.6 would be pending, its
G interval score and cell failures gone, and the look would be unresolved between C0 and C4 at 0.6; the other four
still fail. Those cells were fixed before any fit, and the study was built to put thin priors at risk on them.
