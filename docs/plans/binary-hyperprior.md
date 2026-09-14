# The binary node hyperprior: a default for now

This note answers a ruling of 2026-09-13:

> The reason that motivated the hyperprior was that BART had poor coverage
> on probit models. That could have been poor mixing. Let's extend the
> study, find as reasonable a default as we can for now, and make sure we
> revisit it after we pursue our mixing line of inquiry.

Poor mixing explains almost all of the poor coverage for a sampled k and
almost none of it for a fixed k. The coverage criterion now has an interior
optimum, but a shallow one sitting where chi(1.5, 2) already is. The
earlier real-data result against the hyperprior reverses once that column
is widened. chi(1.5, 2) should stay.

## What was run

Twenty-eight priors: chi(df, scale) for df in {1, 1.25, 1.5, 2, 3} crossed
with scale in {1, 2, 5, Inf}, plus df in {1.5, 3} crossed with scale in
{0.5, 0.25}, and fixed k in {1, 1.5, 2, 3}. Within a case and a repetition
every arm sees the same data and the same seed, so a difference between
arms is a difference in the prior.

Three chain lengths. Short is one chain of 500 draws after 500 discarded -
the earlier study's configuration - on 162 simulated cells and 22 real
datasets. Long is four chains of 2000 after 2000, on 108 simulated cells
(n in {500, 2000}) and 17 real datasets. Probe is eight chains of 8000
after 8000, on four simulated cells and a six-arm subset. Every fit records
held-out log score and Brier, plus split-Rhat and effective sample size for
the sampled k and for the held-out mean probability; the simulated fits
also record the coverage and width of the 90 percent posterior interval for
a held-out row's known true probability. In total 92,160 fits, 65
core-hours.

## Mixing against prior

Take the coverage shortfall of an arm at the short length, 0.90 minus what
it covers, and ask how much of it a longer chain recovers. Over the 540
paired simulated fits per arm:

| arm | short | long | shortfall closed | long coverage |
|---|---|---|---|---|
| chi(1, 1) | 0.828 | 0.930 | 141% | over-covers by 0.030 |
| chi(1.5, 2) | 0.827 | 0.924 | 132% | over-covers by 0.024 |
| chi(1.5, 0.25) | 0.812 | 0.902 | 102% | over-covers by 0.002 |
| k = 1 | 0.792 | 0.855 | 58% | short by 0.045 |
| k = 1.5 | 0.747 | 0.778 | 20% | short by 0.122 |
| k = 2 | 0.679 | 0.700 | 9% | short by 0.200 |
| k = 3 | 0.564 | 0.574 | 3% | short by 0.326 |

That is the hypothesis, confirmed for the hyperprior and refuted for a
fixed k. Every hyperprior arm's shortfall is a chain-length artefact - the
long chains do not merely close it, they overshoot into mild over-coverage.
The fixed arms keep essentially all of theirs.

Going further does little. On the three cells run at all three lengths,
chi(1.5, 2) covers 0.792 short, 0.880 long and 0.892 at the probe: 81
percent of the shortfall closes by the long length and another 11 percent
at the probe, a move of 0.012 with a standard error of 0.006, leaving 0.008
outstanding. chi(1.5, 0.5) moves 0.807 to 0.890 to 0.891 and chi(1.25, 1)
0.814 to 0.888 to 0.890, one and two percent of the shortfall between long
and probe. Fixed k = 2 goes 0.719 to 0.754 to 0.761 and still sits 0.139
short. The hyperprior's coverage has settled by the long length.

The fourth probe cell is the clearest case. In the near-separable cell at
n = 100 with 50 predictors, chi(1.5, 2) covers 0.401 short and 0.956 at the
probe. Fixed k = 2 covers 0.156 at both: a hundred and twenty-eight times
the draws moves it by nothing.

No fixed k reaches nominal coverage at any length. The best at the long
length is k = 1 at 0.855, and its average hides the shape: 0.95 to 0.96 on
the Friedman, linear and weak-signal processes, 0.79 on the strong-signal
one, 0.550 on the near-separable one. Fixed k = 2 reaches 0.90 in 37
percent of the 108 cells against the incumbent's 78 percent.

The sampled k itself, matched on cells and arms:

| length | draws per fit | k split-Rhat | k ESS | fits with k Rhat over 1.05 | probability Rhat | probability ESS |
|---|---|---|---|---|---|---|
| short | 500 | 1.25 | 6 | 69% | 1.034 | 52 |
| long | 8,000 | 1.34 | 12 | 100% | 1.026 | 247 |
| probe | 64,000 | 1.20 | 29 | 98% | 1.016 | 586 |

The sampled k does not converge at any length this study can pay for. A
hundred and twenty-eight-fold increase in draws buys a five-fold increase
in its effective sample size and no convergence; its split-Rhat sits
between 1.2 and 1.4 throughout, and the best single fit at the probe
reached 129. The held-out probability is close to converged at the long
length and fully so at the probe.

Coverage depends on the second and not the first. Across fits it correlates
with k's effective sample size at 0.22, but that is a correlation between
cells rather than within them: remove each cell's mean and it is -0.01 over
540 fits. Cells in which k mixes well are cells that are easy. The coverage
story is about the forest, which the long length fixes, not about k, which
no length fixes.

## Scale below 1

The earlier study reported that within the hyperprior family covering
better was simply being wider - the rank correlation between mean width and
mean coverage across arms was 0.99 - so the coverage criterion named no
interior optimum and would have sent the default off the bottom of the
grid. On converged chains that relation breaks: the same rank correlation
is 0.53 over the 24 hyperprior arms and 0.26 over the 20 finite-scale ones.

| arm | prior median k | sampled k | coverage | coverage distance | distance vs incumbent (SE) | log score vs incumbent (SE) | width |
|---|---|---|---|---|---|---|---|
| chi(1.5, 0.25) | 0.24 | 0.79 | 0.902 | 0.0829 | -0.0026 (0.0036) | +0.0052 (0.0004) | 0.323 |
| chi(1.5, 0.5) | 0.48 | 1.07 | 0.926 | 0.0803 | -0.0052 (0.0022) | +0.0016 (0.0002) | 0.301 |
| chi(1.5, 1) | 0.95 | 1.30 | 0.928 | 0.0815 | -0.0040 (0.0018) | +0.0005 (0.0001) | 0.287 |
| chi(1.5, 2) | 1.91 | 1.48 | 0.924 | 0.0855 | 0 | 0 | 0.279 |
| chi(1.5, 5) | 4.77 | 1.61 | 0.926 | 0.0836 | -0.0019 (0.0016) | -0.0002 (0.0001) | 0.275 |
| chi(3, 0.25) | 0.38 | 0.82 | 0.907 | 0.0819 | -0.0036 (0.0033) | +0.0049 (0.0003) | 0.321 |
| chi(3, 0.5) | 0.77 | 1.11 | 0.926 | 0.0820 | -0.0035 (0.0018) | +0.0013 (0.0002) | 0.299 |
| chi(3, 1) | 1.54 | 1.37 | 0.924 | 0.0855 | -0.0001 (0.0019) | +0.0000 (0.0001) | 0.284 |

Coverage distance is the average over fits of the absolute gap from 0.90.
So there is an interior optimum, at scale 0.5 to 1. Below it the arms get
wider and cover worse: chi(1.5, 0.5) has intervals 0.021 narrower than
chi(1.5, 0.25) (17.8 standard errors) and covers 0.023 more (7.3 standard
errors). The wider-is-better-covering relation inverts.

Four things keep that optimum from deciding anything. It is shallow: the
best improvement anywhere on the grid is 0.0052 of coverage distance on a
base of 0.0855, half a point of coverage on a nominal 90. It is not
resolved at its own bottom - chi(1.5, 0.5) against chi(1.5, 0.25) is 0.9
standard errors. It moves with sample size: at n = 500 scale 0.25 is no
better than the incumbent (+0.0017, SE 0.0055) and at n = 2000 it is the
best arm on the ladder (-0.0092, SE 0.0036). And it is paid for elsewhere:
chi(1.5, 0.5) costs 0.0016 nats of log score (8.4 standard errors) and
doubles the worst-cell log-score regret from 0.0082 to 0.0156, while
chi(1.5, 0.25) costs 0.0052 nats and quadruples it to 0.0354, worse than
fixed k = 2 manages.

Log score has no interior optimum below scale 1 at all: it is flat, 0.2724
to 0.2731, for every arm whose sampled median k stays above about 1.27, and
degrades steadily below - 0.2738 and 0.2742 at a sampled k near 1.1, 0.2774
and 0.2777 near 0.8. Signed coverage at scale 0.25 does land on nominal,
0.9020 against the incumbent's 0.9236, but that is an average of more
over-coverage and more under-coverage rather than better calibration:
15.6 percent of its fits cover below 0.80, against 12.4 percent for the
incumbent.

One structural fact underlies all of this: a thirty-two-fold range of prior
medians, 0.24 to 7.7, produces only a 2.2-fold range in the median k the
sampler draws, 0.79 to 1.73. The data dominate the prior over the whole
grid, which is why every finite-scale arm scores the same.

## Real data, twenty-two datasets

The earlier study's uncomfortable result was that its six real datasets
preferred k held fixed, at about four thousandths of a nat, and it named
this the live question. Sixteen UCI datasets answer it. Held-out log score,
fixed arm minus chi(1.5, 2), so positive is the fixed arm losing:

| arm | 22 datasets, short (SE) | 17 datasets, long (SE) |
|---|---|---|
| k = 1 | +0.0085 (0.0008) | +0.0132 (0.0009) |
| k = 1.5 | +0.0098 (0.0011) | +0.0164 (0.0014) |
| k = 2 | +0.0145 (0.0014) | +0.0226 (0.0019) |
| k = 3 | +0.0268 (0.0022) | +0.0382 (0.0029) |

Brier agrees throughout; fixed k = 2 loses by 0.0039 (SE 0.0004) short and
0.0059 (SE 0.0005) long. By source:

| group | k = 2 against the incumbent, short | long |
|---|---|---|
| the six R datasets | -0.0044 (0.0009) | -0.0027 (0.0007) |
| the sixteen UCI datasets | +0.0216 (0.0019) | +0.0364 (0.0027, eleven) |

The earlier result reproduces exactly on its own six datasets and reverses,
five-fold, on the wider set. Longer chains sharpen it, because the
hyperprior arms gain more from length on real data than the fixed arms do.

What distinguishes the datasets is not size - the correlation between a
dataset's difference and its log training size is 0.00, and the three
largest datasets show differences under 0.001 - and not factors: ten of
the 22 carry factor predictors and their median difference is -0.0004
against +0.0031 for the twelve all-numeric ones. It is separability. Split
the 22 by the median k the incumbent samples: on the ten where it falls
below 1.2 - the well-separated datasets, mean held-out log score 0.15 -
fixed k = 2 loses by 0.034 on average; on the twelve where it exceeds 1.2,
mean log score 0.44, it wins by 0.0018. That is the failure the simulated grid
reports in its near-separable cells - forcing k to 2 over-shrinks a surface
that is nearly deterministic.

The pooled figure is not a typical dataset. On a per-dataset median the
split is eleven to eleven and the median difference is +0.0004; what
decides is the asymmetry of the tails. The hyperprior's worst dataset costs
it 0.0091 against fixed k = 2 (kyphosis), while fixed k = 2's worst costs
it 0.184 (tic-tac-toe), 0.079 (sonar), 0.051 (ionosphere) and 0.033
(banknote) at the long length. That twenty-fold asymmetry is the argument,
and it is the worst-case argument that set the default originally.

## The default for now

The decision rule, stated before it is applied. Move the shipped default
only if some arm improves the criterion that motivated the hyperprior -
coverage of the 90 percent interval, on converged chains - by an amount a
user would see in a reported coverage figure, which is two points, or 0.02
in mean coverage distance; and only if it does not lose noticeably on
anything else, taking a hundredth of a nat as noticeable on log score,
simulated or real, and a doubling as noticeable on worst-cell regret. If
nothing clears the coverage bar, keep chi(1.5, 2): a default with a history
should not move for a difference the next study would not reproduce.

| arm | coverage | distance vs incumbent (SE) | simulated log score (SE) | worst-cell regret | real log score, 22 (SE) |
|---|---|---|---|---|---|
| chi(1.5, 2) | 0.924 | 0 | 0 | 0.0082 | 0 |
| chi(1.5, 5) | 0.925 | -0.0019 (0.0016) | -0.0002 (0.0001) | 0.0071 | +0.0004 (0.0004) |
| chi(1.5, 1) | 0.928 | -0.0040 (0.0018) | +0.0005 (0.0001) | 0.0093 | -0.0003 (0.0005) |
| chi(1, 1) | 0.930 | -0.0052 (0.0020) | +0.0005 (0.0001) | 0.0127 | +0.0003 (0.0005) |
| chi(1.5, 0.5) | 0.925 | -0.0052 (0.0022) | +0.0016 (0.0002) | 0.0156 | +0.0021 (0.0006) |
| chi(1.5, 0.25) | 0.902 | -0.0026 (0.0036) | +0.0052 (0.0004) | 0.0354 | +0.0068 (0.0009) |
| k = 1 | 0.855 | +0.0418 (0.0049) | +0.0041 (0.0004) | 0.0390 | +0.0085 (0.0008) |
| k = 2 | 0.700 | +0.1666 (0.0096) | +0.0054 (0.0004) | 0.0448 | +0.0145 (0.0014) |

Nothing clears the bar. The largest coverage improvement available anywhere
in the grid is 0.0052, a quarter of the threshold, and the two arms that
offer it both pay for it - chi(1, 1) with a 55 percent rise in worst-cell
regret, chi(1.5, 0.5) with a doubling of it plus 0.0021 nats on real data.
Keep chi(1.5, 2).

What has changed is the reason. The earlier note kept the incumbent because
the coverage criterion had no stopping point; this study finds the stopping
point, and it is where the incumbent already sits, to within half a point
of coverage. The grid is flat on every point-prediction score for a prior
median between about 0.8 and 7.7 and degrades below that, and
chi(1.5, 2)'s prior median of 1.91 is comfortably inside. Both earlier
negatives survive and one strengthens: the improper scale is worse at every
degree of freedom, by 0.012 to 0.026 of coverage distance at df 1.5 and
above, and no fixed k is tolerable, now on real data as well as simulated.

## What to revisit after the mixing work

Three results, and only three, rest on the sampled k's non-convergence.
First, every sampled-k number quoted here is a chain-length-dependent
summary rather than a posterior one; they are comparable across arms
because every arm ran identically, and that is all the weight they carry.
The appendix's per-dataset k medians in particular are labels for
separability, not estimates. Second, the residual 0.008 of coverage
shortfall the hyperprior still shows at the probe, and the process-level
residues behind it - chi(1.5, 2) reaches only 0.873 on the strong-signal
cells and 0.862 on the near-separable ones at the long length - could be
the sampler or the model. Third, and most important for the default, the
interior optimum is located by coverage-distance differences of 0.002 to
0.005 between adjacent scales, the same order as the 0.012 the incumbent's
own coverage still moved between the long length and the probe. The
optimum's location is not resolved to better than about a factor of two in
scale.

To reopen the default, a mixing fix would have to take the sampled k to a
split-Rhat below 1.05 with an effective sample size in the hundreds at a
length a user would actually run - against the probe's 64,000 draws per fit
reaching 29 - and then, at that length, produce either a spread in coverage
distance across the finite-scale arms larger than 0.02 or an optimum
somewhere other than scale 0.5 to 2. If the spread stays at half a point of
coverage, no mixing fix reopens the question; it only makes this verdict
better founded. The fixed-k results, the largest effects in this note,
involve no sampled k at all.

## What this study cannot settle

- Whether the coverage residue in the strong-signal and near-separable
  cells is the leaf prior, the tree count, the probit link or the forest.
  It survives a 128-fold increase in draws, so it is not mixing, and
  nothing here separates the rest.
- Coverage at n = 100 on converged chains. The long leg drops that size and
  only one n = 100 cell appears in the probe.
- The sampled k as a posterior quantity, at any length.
- Large real datasets at length. Five of the 22 exceed 2,000 rows and have
  only the short leg, and every dataset over 5,000 rows is subsampled to
  4,000 training rows per split.
- The ordering among finite-scale arms on any score: within noise
  everywhere, and it changes from score to score.
- Anything off the fitted configuration of 75 trees and package defaults.
  The tree count interacts with k by construction.
- The logistic and weighted binary paths. Only probit was fitted.

## Re-running

    Rscript benchmarks/R/binary-hyperprior.R blocks
    Rscript benchmarks/R/binary-hyperprior.R sim:friedman:500 <outdir>
    Rscript benchmarks/R/binary-hyperprior.R sim:strong:2000:20:0.2 <outdir>
    Rscript benchmarks/R/binary-hyperprior.R real:sonar <outdir>
    Rscript benchmarks/R/binary-hyperprior.R summarize <outdir>
    Rscript benchmarks/R/binary-hyperprior.R compare <shortdir> <longdir>

One block is one invocation and writes one rds, so a run splits across a
session and a restart costs only the interrupted block. Appending a
predictor count or a base rate narrows a simulated block, which is how the
long and probe legs were kept affordable. `compare` pairs two directories
on cell, repetition and arm, so with a short directory first it reports
what the chain length does.

The legs differ only in environment variables:

| leg | variables |
|---|---|
| short | `BINARY_HYPERPRIOR_CHAINS=1 _BURN=500 _DRAWS=500` |
| long | `BINARY_HYPERPRIOR_CHAINS=4 _BURN=2000 _DRAWS=2000` |
| probe | `BINARY_HYPERPRIOR_CHAINS=8 _BURN=8000 _DRAWS=8000 BINARY_HYPERPRIOR_ARMS=mixing` |

`BINARY_HYPERPRIOR_REPS` sets the simulated repetitions per cell (8 at
n = 100 and n = 500 short, 6 at n = 2000 short and n = 500 long, 4 at
n = 2000 long and at the probe), `BINARY_HYPERPRIOR_SPLITS` the real-data
splits (40 throughout), `BINARY_HYPERPRIOR_ARMS` an arm list or the named
`mixing` subset, and `BINARY_HYPERPRIOR_CORES` the worker count. The
settings ride each rds, so `summarize` reports which length produced a
directory and refuses to mix incompatible blocks. UCI files download on
first use into `DBARTS_BENCH_DATA`, or the package's user cache when that
is unset, and are checked against a recorded sha256.

Cost in core-hours: five for the short leg, two for the probe, 44 for the
simulated long leg, 14 for the UCI legs.

## Appendix: the twenty-two real datasets

Rows and predictors are after cleaning; "train" is the training rows per
80/20 split, capped at 4,000; "factors" is how many predictors are factors.
"k" and "log score" are the median k chi(1.5, 2) samples and the log score
it attains; k orders the table because it is what predicts the result. The
last column is fixed k = 2 minus chi(1.5, 2) on held-out log score at the
short length, positive meaning fixed k loses. An asterisk marks the
datasets that also ran at the long length.

| dataset | source | rows | pred | rate | train | factors | log score | k | k = 2 - chi |
|---|---|---|---|---|---|---|---|---|---|
| tictactoe* | UCI | 958 | 9 | 0.653 | 766 | 9 | 0.066 | 0.32 | +0.1845 |
| banknote* | UCI | 1372 | 4 | 0.445 | 1097 | - | 0.018 | 0.57 | +0.0283 |
| mushroom | UCI | 8124 | 21 | 0.482 | 4000 | 21 | 0.003 | 0.59 | +0.0060 |
| spambase | UCI | 4601 | 57 | 0.394 | 3680 | - | 0.193 | 0.63 | +0.0037 |
| ionosphere* | UCI | 351 | 33 | 0.641 | 280 | - | 0.185 | 0.68 | +0.0391 |
| sonar* | UCI | 208 | 60 | 0.534 | 166 | - | 0.373 | 0.76 | +0.0622 |
| wdbc* | UCI | 569 | 30 | 0.373 | 455 | - | 0.101 | 0.82 | +0.0142 |
| climate* | UCI | 540 | 18 | 0.085 | 432 | - | 0.148 | 0.84 | +0.0076 |
| biopsy* | R | 683 | 9 | 0.350 | 546 | - | 0.096 | 1.04 | -0.0040 |
| magic | UCI | 19020 | 10 | 0.352 | 4000 | - | 0.343 | 1.12 | -0.0007 |
| bank | UCI | 45211 | 16 | 0.117 | 4000 | 9 | 0.224 | 1.21 | -0.0005 |
| adult | UCI | 48842 | 14 | 0.239 | 4000 | 8 | 0.303 | 1.22 | -0.0004 |
| creditapproval* | UCI | 653 | 15 | 0.453 | 522 | 9 | 0.329 | 1.45 | +0.0018 |
| cleveland* | UCI | 297 | 13 | 0.461 | 237 | 4 | 0.364 | 1.68 | -0.0029 |
| pbc* | R | 276 | 16 | 0.402 | 220 | 1 | 0.490 | 1.69 | -0.0066 |
| kyphosis* | R | 81 | 3 | 0.210 | 64 | - | 0.397 | 1.86 | -0.0091 |
| pima* | R | 532 | 7 | 0.333 | 425 | - | 0.455 | 1.87 | -0.0025 |
| german* | UCI | 1000 | 20 | 0.300 | 800 | 13 | 0.485 | 2.16 | -0.0003 |
| transfusion* | UCI | 748 | 4 | 0.238 | 598 | - | 0.478 | 2.53 | +0.0009 |
| infert* | R | 248 | 5 | 0.335 | 198 | 1 | 0.589 | 2.73 | -0.0044 |
| birthwt* | R | 189 | 8 | 0.312 | 151 | 2 | 0.588 | 2.78 | -0.0001 |
| haberman* | UCI | 306 | 3 | 0.265 | 244 | - | 0.546 | 2.82 | +0.0025 |
