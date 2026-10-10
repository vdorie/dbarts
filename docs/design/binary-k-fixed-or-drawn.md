# Binary fits: should k be fixed or drawn?

Status: FINDING, 2026-10-09 (dec-B404, dec-B390); ruled dec-B420. The September binary
hyperprior study rerun with the probit rescaling step in place, on probit
and logistic, fixed k against drawn k. The per-fit results are in
scratch/khp-study, untracked.

## Summary

k is the leaf prior's scale. Binary fits draw it by default under chi(1.5,
2). The question was whether drawing it is worth anything, or whether a
fixed value would do as well.

**A fixed k is not good enough, under probit or logistic.** The best fixed
value, k = 1, misses four of the six bars set below in both links. It
covers the true probability worse, it loses on real data and it has a
worst case eight times the default's. The loss comes from two directions at
once. On nearly separated data k = 1 is too large: tictactoe loses 0.09 to
0.10 nats. On weak-signal data it is too small: haberman loses 0.03 and
birthwt 0.02, where the drawn k sits near 2.8. No single fixed value suits
both, and that is what the drawn k buys.

**chi(1.5, 2) stays the default for both links.** No drawn prior improves
coverage by the two points the September rule requires, at either chain
length. The closest, chi(1.5, 0.5) under logistic, gains 1.65 points of
coverage distance. It pays 0.003 nats on real data and more than doubles
the worst-case regret.

**The logistic collapsed scale move should be planned.** dec-B404 builds it
if a drawn k is kept, and a drawn k is kept. Under logistic, k does not
converge at any length tried: at 64,000 draws every fit still has a split-Rhat
above 1.05, and k is still drifting. What a user sees is mostly unaffected,
since the predicted probabilities and log score are converged. Two things do
depend on the chain length: the reported k, and the intervals for
probabilities very near 0 or 1.

One finding bears on probit too. At 75 trees, n = 500 and up to 50
predictors, the probit rescaling step does not converge k either. Even at
64,000 draws k's split-Rhat is 1.12 to 1.21. September's condition for
reopening the default was a split-Rhat below 1.05 with an effective sample
size in the hundreds. The step does not meet that condition at these sizes.

## What was run

The September study's harness, benchmarks/R/binary-hyperprior.R, with
three additions: a setting that fits the logistic link to the same cases,
an arm subset for this study, and a second coverage figure, described
under "Interior coverage" below. 75 trees and package defaults throughout.

Nine arms:

- k fixed at 1, 1.5, 2 and 3;
- the default chi(1.5, 2);
- the four drawn priors closest to it in September's table: chi(1.5, 5),
  chi(1.5, 1), chi(1, 1) and chi(1.5, 0.5).

Within a case every arm and both links see the same data and the same
seed. The simulated truth is a probit surface for both links.

Two chain lengths:

- **Long** is 4 chains of 2000 draws after 2000. It ran on all 108
  simulated cells at n = 500 and 2000: 6 processes, 3 predictor counts and 3
  base rates. It also ran on the 22 real datasets.
- **Probe** is 8 chains of 8000 draws after 8000, for the five drawn arms
  only. It ran on the 9 near-separable cells at n = 500 and the 3
  strong-signal cells at n = 500 with 20 predictors. These are where
  September found the differences. The fixed arms are compared at the long
  length, since their predictions are converged there.

Two additions:

- the near-separable cells at n = 100, at the long length;
- the drawn arms on tictactoe at the probe length, 10 splits.

The full design came to about 73 core-hours against a 48 core-hour budget,
so repetitions were cut:

- simulated cells: 4 repetitions at n = 500 and 3 at n = 2000, against
  September's 6 and 4. They use September's first seeds.
- real datasets: 20 splits, or 10 on the five with 3,700 to 4,000 training
  rows (mushroom, spambase, magic, bank, adult), against September's 40.

In all, 13,824 long fits, 480 probe fits and 648 at n = 100.

Per fit at the long length, on one x86-64 core, the times were as follows.

| rows | probit | logistic |
|---|---|---|
| 500 | 6.2 s | 8.0 s |
| 2000 | 9.3 s | 16.7 s |
| 4000 | 19 s | 33 s |

The arms within a link always ran on one machine, so each pairing is exact.
The probit real data and the two additions ran on an arm64 Mac; the rest
ran on the x86-64 host.

## The bars

The decision rule, stated before any result was seen:

- Coverage distance is the average over fits of |coverage of the 90 percent
  interval - 0.90|. A difference is noticeable at 0.02. This is
  September's bar, from docs/plans/binary-hyperprior.md.
- Log score, simulated or real: noticeable at 0.01 nats. This is
  September's bar.
- Worst-cell regret on log score: noticeable at a doubling. This is
  September's bar.
- Brier, simulated or real: noticeable at 0.0025. September set no Brier
  bar. This one was set for this study, at a quarter of the log-score bar,
  the ratio September's real-data differences showed.

A fixed k is good enough if it is within every bar of the best drawn prior
on every criterion. The default moves only if some arm improves coverage
distance by 0.02 without a noticeable loss elsewhere.

## Fixed against drawn

Long length. Each entry is k fixed minus the best drawn arm on that
criterion; the last column is the fixed arm's worst-cell regret over the
best drawn arm's.

| link | arm | coverage distance | sim log score | sim Brier | real log score | real Brier | worst-cell regret |
|---|---|---|---|---|---|---|---|
| probit | k = 1 | +0.050 | +0.0034 | +0.0009 | +0.0124 (0.0014) | +0.0026 | 8.0 times |
| probit | k = 1.5 | +0.116 | +0.0024 | +0.0004 | +0.0158 | +0.0035 | 12.5 times |
| probit | k = 2 | +0.174 | +0.0044 | +0.0008 | +0.0217 | +0.0054 | 16.6 times |
| probit | k = 3 | +0.282 | +0.0113 | +0.0024 | +0.0362 | +0.0102 | 24.4 times |
| logistic | k = 1 | +0.070 | +0.0041 | +0.0012 | +0.0125 (0.0012) | +0.0028 | 7.9 times |
| logistic | k = 1.5 | +0.129 | +0.0036 | +0.0007 | +0.0155 | +0.0035 | 12.5 times |
| logistic | k = 2 | +0.185 | +0.0056 | +0.0010 | +0.0210 | +0.0051 | 16.9 times |
| logistic | k = 3 | +0.270 | +0.0125 | +0.0025 | +0.0346 | +0.0095 | 24.9 times |

Standard errors are paired, in parentheses where they bear on a bar. Every
fixed arm's worst cell is tictactoe in both links.

k = 1 misses coverage distance and worst-cell regret by a wide margin. It
misses real log score and Brier narrowly.

On real data, the per-dataset log score of k = 1 minus the default:

| | probit | logistic |
|---|---|---|
| all 22 datasets | +0.0124 | +0.0125 |
| without tictactoe | +0.0077 | +0.0084 |
| median dataset | +0.0062 | +0.0063 |
| datasets losing more than 0.01 | 6 of 22 | 6 of 22 |
| tictactoe | +0.100 | +0.088 |
| haberman | +0.030 | +0.032 |
| sonar | +0.024 | +0.028 |
| birthwt | +0.016 | +0.019 |
| banknote | +0.017 | +0.015 |
| ionosphere | +0.013 | +0.014 |

A longer chain widens the gap. At the probe length the drawn arms gain
about 0.005 nats on tictactoe and on the near-separable cells, while the
fixed arms are already converged at the long length.

| link | drawn arms on tictactoe, long | drawn arms on tictactoe, probe | k = 1 |
|---|---|---|---|
| probit | 0.031 to 0.036 | 0.028 to 0.030 | 0.133 |
| logistic | 0.041 to 0.043 | 0.036 | 0.131 |

## The default

Long length. Differences are from chi(1.5, 2), with paired standard errors.

| link | arm | coverage | coverage distance | sim log score | real log score | worst-cell regret |
|---|---|---|---|---|---|---|
| probit | chi(1.5, 2) | 0.928 | 0 | 0 | 0 | 0.0155 |
| probit | chi(1.5, 5) | 0.923 | +0.0039 (0.0021) | -0.0001 | +0.0010 | 0.0128 |
| probit | chi(1.5, 1) | 0.930 | -0.0004 (0.0026) | +0.0006 | +0.0002 | 0.0138 |
| probit | chi(1, 1) | 0.913 | +0.0107 (0.0032) | +0.0016 | +0.0024 | 0.0207 |
| probit | chi(1.5, 0.5) | 0.914 | +0.0026 (0.0032) | +0.0025 | +0.0030 | 0.0270 |
| logistic | chi(1.5, 2) | 0.794 | 0 | 0 | 0 | 0.0113 |
| logistic | chi(1.5, 5) | 0.791 | +0.0018 (0.0021) | -0.0001 | +0.0008 | 0.0153 |
| logistic | chi(1.5, 1) | 0.798 | -0.0037 (0.0022) | +0.0003 | +0.0008 | 0.0146 |
| logistic | chi(1, 1) | 0.809 | -0.0105 (0.0024) | +0.0003 | +0.0007 | 0.0162 |
| logistic | chi(1.5, 0.5) | 0.816 | -0.0165 (0.0027) | +0.0015 | +0.0029 | 0.0254 |

Under probit nothing comes near the 0.02 coverage bar. Under logistic
chi(1.5, 0.5) comes closest. It stays under the bar and pays 0.0029 nats
on real data and 2.2 times the worst-cell regret.

At the probe length, on its 12 cells, nothing clears the bar either. The
largest improvement is chi(1, 1) under logistic, -0.021 (0.013) on 48
fits per arm. On interior coverage the same arm is -0.0007.

## Interior coverage

Logistic covers 0.79 on the simulated cells where probit covers 0.93, and
0.53 against 0.89 on the near-separable ones. That gap is not the prior.

The simulated truth is a probit surface. On the near-separable cells only
20 percent of held-out rows have a true probability between 0.01 and 0.99.
Past that range the probit truth goes to 0 or 1 far faster than a logistic
fit can follow, so a logistic interval misses the true probability while
being off by a tiny absolute amount.

A second coverage figure was therefore recorded: coverage over the rows
whose true probability is within [0.01, 0.99]. It was added after a first
look at the data, so it is a sensitivity reading, not part of the rule.

| link | arm | interior coverage, all 108 cells | interior coverage, near-separable | full coverage, near-separable |
|---|---|---|---|---|
| probit | chi(1.5, 2) | 0.951 | 0.956 | 0.893 |
| probit | chi(1.5, 0.5) | 0.941 | 0.942 | 0.846 |
| probit | k = 1 | 0.955 | 0.980 | 0.551 |
| probit | k = 2 | 0.913 | 0.831 | 0.204 |
| logistic | chi(1.5, 2) | 0.956 | 0.977 | 0.525 |
| logistic | chi(1.5, 0.5) | 0.951 | 0.975 | 0.546 |
| logistic | k = 1 | 0.955 | 0.980 | 0.277 |
| logistic | k = 2 | 0.907 | 0.823 | 0.168 |

Read this way, every drawn arm is within 0.0014 of the default's interior
coverage distance in both links. Fixed k = 1 is as good as any drawn arm.
k = 1's whole coverage shortfall lies in the near-0 and near-1 rows of the
near-separable and strong-signal cells.

At n = 100 on the near-separable cells, where only 20 percent of rows are
interior, k = 1 trails on interior coverage too:

| link | chi(1.5, 2) | k = 1 |
|---|---|---|
| probit | 0.914 | 0.883 |
| logistic | 0.907 | 0.881 |

The verdict on a fixed k does not rest on full coverage. With interior
coverage in its place, k = 1 still misses on real log score, real Brier and
worst-cell regret.

## Convergence

Long length:

| link | arms | k split-Rhat above 1.05 | k effective draws, median (10th percentile) | probability split-Rhat above 1.05 | probability effective draws, median (10th percentile) |
|---|---|---|---|---|---|
| probit | drawn | 68 to 72% of fits | 26 to 28 (7) | 30 to 36% | 235 to 384 (10 to 14) |
| probit | k = 1 | - | - | 7% | 376 (101) |
| logistic | drawn | 90 to 93% | 12 to 13 (6) | 6 to 7% | 578 to 756 (99 to 110) |
| logistic | k = 1 | - | - | 3% | 883 (227) |

"Probability" here is the held-out mean probability, the forest's own
summary.

Probe against long, drawn arms, ranges over the five:

| link, cells | coverage, long to probe | interior coverage, long to probe | k split-Rhat at probe | k effective draws at probe | probability split-Rhat above 1.05 at probe |
|---|---|---|---|---|---|
| probit, near-separable | 0.77-0.89 to 0.81-0.91 | 0.93-0.96 to 0.96-0.98 | 1.13 to 1.21 | 38 to 57 | 69 to 78% |
| probit, strong | 0.88-0.92 to 0.90-0.93 | 0.90-0.93 to 0.92-0.94 | 1.12 to 1.27 | 29 to 63 | 25 to 67% |
| logistic, near-separable | 0.57-0.61 to 0.72-0.77 | 0.98-0.99, unchanged | 1.55 to 1.64 | 15 to 16 | none |
| logistic, strong | 0.72-0.76 to 0.80-0.86 | 0.95, unchanged | 1.41 to 1.62 | 15 to 21 | none |

Under logistic at the probe length:

- k's split-Rhat exceeds 1.05 in every fit;
- k's median on the near-separable cells moves from 0.44 to 0.35 between
  the two lengths;
- on tictactoe it moves from 0.16 to 0.11.

The predicted probabilities are converged (split-Rhat 1.01, 1,300 to 2,400
effective draws), and log score moves by 0.001 with length. Full coverage
rises by 0.08 to 0.16 with length, because the near-0 and near-1
probabilities follow k.

Under probit the rescaling step improves on September's sampler: k's
effective draws at 64,000 draws rise from about 29 to 38 to 63. It still
leaves both k and, on the near-separable cells, the predicted probability
unconverged. The step's own measurements (probit-k-scale-move.md) used 150
rows, 3 predictors and 50 trees.

What is a posterior quantity, and what is a summary of the chain length:

- **Chain-length summaries:** every sampled-k figure, in both links and at
  both lengths. Probit's drawn-arm coverage on near-separable and
  well-separated cases. Logistic's coverage of near-0 and near-1 rows.
- **Posterior quantities:** every fixed-arm figure. Logistic's predicted
  probabilities, log score, Brier and interior coverage. Probit's log score
  and Brier, which move by 0.005 or less between lengths.

The comparisons between arms are fair either way, since every arm ran
identically, and the longer chain only widened the drawn arms' lead.

## For the logistic move

k-mixing-pg-families.md found logistic's k slow at small true k and near
separation. k-tail-samplers.md measured the collapsed scale move taking
logistic's 99th-percentile autocorrelation time from about 23,000 sweeps
to 66. This study adds three things:

- a drawn k is worth keeping under logistic, so the move has a k to serve;
- at the sizes users fit, logistic's k is unconverged at 8,000 and 64,000
  draws;
- the user-visible cost of that is the reported k and the intervals near 0
  and 1, not the predictions.

## What this cannot settle

- **The Brier bar.** It was set for this study, not by the maintainer. k = 1
  misses it by 0.0001 (probit) and 0.0003 (logistic).
- **Dependence on tictactoe.** k = 1's real-data log-score miss depends on
  tictactoe. Without it the loss is 0.0077 and 0.0084, under the bar. The
  worst-cell regret and the two-sided losses (haberman, birthwt) do not
  depend on it.
- **k as a posterior quantity.** It is not settled under either link at any
  length tried. The ordering among the drawn priors rests on chains whose k
  has not converged. It agrees with September's and with this study's probe.
- **Probit at these sizes.** The probit rescaling step does not meet
  September's reopen condition at n = 500, 75 trees and up to 50
  predictors. Whether the rest needs the collapsed move as well
  (k-tail-samplers.md) or a tree-level move (TODO tree-mixing-proposals)
  is not settled here.
- **Logistic on a logistic truth.** Every simulated truth is a probit
  surface, so logistic's full coverage on strong signals measures the two
  links' tails as much as the fit.
- **Large real datasets.** The five datasets of 3,700 to 4,000 training
  rows ran 10 splits.
- **Other configurations.** Anything off 75 trees and package defaults.

Earlier notes: probit-k-calibration.md (k's slow mixing under probit),
probit-k-scale-move.md (the probit rescaling step), k-mixing-pg-families.md
(k under logistic and negative binomial), k-tail-samplers.md (samplers for
the remaining tail).
