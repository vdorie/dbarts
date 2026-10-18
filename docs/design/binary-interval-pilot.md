# Binary fits: does the collapsed scale move make intervals better calibrated?

Status: FINDING, 2026-10-17 (dec-B462). The pilot of the binary interval study, stopped there. It compared
the shipped sampler with the same sampler plus the collapsed scale move, under probit and logistic, at equal
running time. The plan, [binary-interval-study.md](../plans/binary-interval-study.md), is closed; a different study
replaces its main run ([binary-k-mixing-need.md](../plans/binary-k-mixing-need.md)). Per-fit results are in
scratch/bis-run, untracked.

## Summary

You asked for proof that the collapsed move helps what users report (intervals, effects, false declarations),
not k itself. The pilot did not show it, and it did not show harm either.

- **Effects and false declarations: no gain worth the bar.** Under probit the effect group moved by +0.0035
  (bounds -0.0025 to +0.0089); under logistic by +0.0085 (+0.0056 to +0.0140). The bar was +0.02. Nothing on
  false declarations or wrong signs moved by more than half a point, except the four-arm design, whose four
  cells are too few to read.
- **Probabilities: the move looks harmful as the rule was written, and that reading is not valid.** The loss
  (2.0 points under probit, 5.9 under logistic) sits in datasets where the true probability is within 0.01 of
  0 or 1 on most rows. There a posterior mean falls outside its own interval, so a sampler that reaches the
  posterior "misses" by construction, and the move has the lower interval score on all three of the worst datasets and the lower error on two. With those
  datasets set apart, the gain is +0.0007 under probit and +0.0019 under logistic, nothing either way.
- **Real outcomes: the apparent gap was mostly an artifact.** Once saturated rows and single-precision
  rounding were handled, default-length runs of both samplers reproduce the long-run intervals about equally
  well (under logistic 0.872 shipped against 0.881 move, nominal 0.90).
- **The move costs 0.6 to 30 percent more time per iteration, 9 to 10 percent on average**, less than the
  14 to 28 percent expected.
- **The pilot cannot say whether slow k mixing hurts the estimands users report.** It was stopped for that
  reason and a different study asks it directly.

## What was compared

Two builds of the same code. The shipped sampler is what bartcore ships: probit with its rescaling step,
logistic with no scale move. The move build adds the collapsed scale move (to logistic; to probit beside the
rescaling step). Package defaults otherwise (75 trees, k drawn under chi(1.5, 2), 500 burn-in). Both builds
were checked on the x86 box before any study fit: the shipped build draws bit for bit what the reference build
draws, on every path the study uses; the move build differs on every binary fit; a longer run's first draws
equal a shorter run's.

**Equal running time.** The move takes longer per iteration. Each cell's time ratio was measured on its own
(below), and the shipped sampler was then given that many more iterations, spent either as more draws or as
more burn-in. Pooled over the 166 cells of the probability, effect and zero-coverage groups (83 a link), the
distance from nominal was 0.1263 with more draws and 0.1280 with more burn-in, so more draws was kept. The
two give nearly the same answers. The move is read at the default length (500 draws after 500 burn-in).

**Truths known.** Real covariates, outcomes drawn from a known truth. 79 cells, 10 replications each, so 790
paired fits a link (a pair is the shipped and the move fit on the same data and seed); no pair was skipped and
no draw was non-finite.

- Probabilities: the 22 real binary datasets of the k study, each with two truths fitted to its real outcome,
  one by BART at long length and one by a smooth additive model (44 cells). Held-out rows are scored against
  the true probability.
- Effects of a treatment on a binary outcome, each unit's effect scored against its true value (35 cells):
  IHDP's covariates, the ACIC 2016 competition's covariates and generated surfaces (including 250 covariates,
  effects in a subgroup only, a 1 percent outcome and weak overlap), the ACIC 2019 binary practice sets,
  Twins, and a four-arm design that stands for the multiple-comparisons case. 27 of these cells have real
  effects; the rest have none.

**Real outcomes, truth unknown.** The 22 datasets, adult at full size and three causal datasets (lindner,
ACTG 175, lalonde) with their real outcomes, 26 problems a link. Each was compared with a reference of 16 long
chains, 8 from each build. Neither build has a known truth here, so this asks only whether a default-length
run reproduces the long-run posterior.

**How a cell is read.** A cell's gap is how far the coverage of its 90 percent intervals is from 0.90, on
either side, corrected for the noise that makes a measured distance look larger than the true one. A group's
gap is the mean over its cells. The move's gain is the shipped gap minus the move's gap: positive means the
move is closer to nominal. 0.02 is 2 points of coverage. The move was to be kept if a group gained at least
0.02 with its lower bound above zero and no guard showed harm; guards are listed beside the readings. Bounds
below are one-sided at 0.005. The pilot spends none of the study's error budget.

## The pilot's readings, frozen rule

Gain of the move over the shipped sampler at equal time; bounds in brackets. Shipped and move columns are
the groups' gaps.

| group (cells, pairs) | link | shipped | move | gain [bounds] |
|---|---|---|---|---|
| Probabilities (44, 440) | probit | 0.1260 | 0.1461 | -0.0200 [-0.0268, -0.0137] |
| | logistic | 0.0658 | 0.1252 | -0.0594 [-0.0676, -0.0516] |
| Effects, real effects (27, 270) | probit | 0.1938 | 0.1904 | +0.0035 [-0.0025, +0.0089] |
| | logistic | 0.2271 | 0.2186 | +0.0085 [+0.0056, +0.0140] |
| Zero effect, coverage of zero (12, 120) | probit | 0.0475 | 0.0450 | +0.0025 [-0.0041, +0.0026] |
| | logistic | 0.0494 | 0.0500 | -0.0007 [-0.0027, 0.0000] |
| Zero effect, any unit declared (11, 110) | probit | 0.0886 | 0.0932 | -0.0045 [-0.0341, +0.0455] |
| | logistic | 0.0500 | 0.0500 | 0 |
| Declared effects with the wrong sign (26, 260) | probit | 0.0304 | 0.0308 | -0.0004 [-0.0041, +0.0033] |
| | logistic | 0.0409 | 0.0368 | +0.0041 [-0.0009, +0.0087] |
| One band over all units (78, 780) | probit | 0.4630 | 0.4831 | -0.0200 [-0.0412, -0.0029] |
| | logistic | 0.4401 | 0.4421 | -0.0021 [-0.0228, +0.0235] |
| Average effect (35, 350) | probit | 0.0898 | 0.0929 | -0.0031 [-0.0158, +0.0134] |
| | logistic | 0.0771 | 0.0796 | -0.0025 [-0.0178, +0.0059] |

Two ratios, move over shipped, upper bound at the same level: root-mean-square error of the posterior mean
1.0009 (1.0135) probit, 1.0026 (1.0105) logistic; interval score 0.9953 (1.0078) probit, 0.9628 (0.9707)
logistic, over all 79 cells. The move never raises the error, and under logistic its interval score is about 4
percent better. The four-arm design alone (4 cells, 40 pairs) has bounds too wide to read: its
any-unit-declared gain under probit is -0.0125 with bounds -0.0937 to +0.1250.

The guards were the largest worsening the move may show before it counts as harm: 0.010 for probabilities,
0.015 for effects and for coverage of zero, 0.075 for any unit declared, 0.010 for the wrong-sign share,
0.035 for the band, 0.040 for the average effect, and 1.02 for each ratio. Every margin stayed at its planned
value after the pilot. Read as written, the rule says "leave out: harm shown" under both links, on the
probability group alone: its lower bound on the worsening (0.0137 probit, 0.0516 logistic) is past 0.010.

## Why the probability reading shows harm

The harm sits in near-deterministic cells: a dataset and truth whose true probability lies inside
[0.01, 0.99] on fewer than half of its rows. The truth is the posterior mean of a fitted probability surface,
and on near-separable data that posterior is so skewed that its mean falls outside its own 90 percent
interval. A sampler that reaches the posterior then fails to cover the truth, however correct it is.

- **The move fits the posterior better there.** Logistic, BART-fitted truth, shipped (equal time) against
  move, coverage of the true probability, interval score (width plus a penalty for misses; lower is better)
  and error:

| dataset | rows with truth in [0.01, 0.99] | coverage | interval score | error |
|---|---|---|---|---|
| mushroom | 0.04 percent | 0.973 vs 0.076 | 0.0092 vs 0.0054 | 0.0079 vs 0.0072 |
| banknote | 2.0 percent | 0.986 vs 0.194 | 0.0633 vs 0.0498 | 0.0341 vs 0.0349 |
| tictactoe | 5.7 percent | 0.919 vs 0.497 | 0.2214 vs 0.1396 | 0.1054 vs 0.0947 |

  On tictactoe's rows with truth inside [0.01, 0.99] the move covers 0.90 against the shipped sampler's 0.59.
- **The mechanism, from the review of these fits.** In the move's own logistic fits to the real outcomes of
  mushroom, banknote and tictactoe, the posterior mean lies outside its own 90 percent interval on 99.3, 95.2
  and 92.7 percent of rows, split about evenly between rows near 0 and near 1. The shipped sampler's share
  rises from 0.235 at the default length to 0.55 at its long reading: it covers when short because it has not
  reached that posterior. The move's interval ends are exactly 0 or 1 in double precision on 47 to 99.8
  percent of rows of these datasets.
- **The additive truths have the same problem.** For banknote, wdbc, climate, ionosphere and sonar not one
  row of the truth lies inside [0.01, 0.99] (mushroom: 0.04 percent), and the truth is held at a floor of
  0.0001 where the data cannot tell 0.0001 from 0.000000000001. Coverage there is set by the floor.
- **What is real.** Under separation a converged drawn-k posterior is extremely concentrated near 0 or 1, and
  the shipped sampler's intervals there depend on chain length. That is the finding; it is not a calibration
  loss.

Classifying each cell and link from its truth alone, 16 cells under probit and 17 under logistic are
near-deterministic (of 44): biopsy, mushroom, banknote, tictactoe, wdbc, climate and ionosphere, both truths;
sonar and spambase on the additive truth; spambase on the BART truth under logistic only.

**Frozen and amended side by side**, with those cells left out of the probability group and of the band's
probability cells:

| reading | link | cells | shipped | move | gain [bounds] | verdict |
|---|---|---|---|---|---|---|
| probabilities, frozen | probit | 44 | 0.1260 | 0.1461 | -0.0200 [-0.0268, -0.0137] | leave out, harm shown |
| probabilities, amended | probit | 28 | 0.0415 | 0.0407 | +0.0007 [-0.0034, +0.0050] | inconclusive |
| band, frozen | probit | 78 | 0.4630 | 0.4831 | -0.0200 [-0.0412, -0.0029] | |
| band, amended | probit | 62 | 0.4198 | 0.4271 | -0.0074 [-0.0275, +0.0096] | |
| probabilities, frozen | logistic | 44 | 0.0658 | 0.1252 | -0.0594 [-0.0676, -0.0516] | leave out, harm shown |
| probabilities, amended | logistic | 27 | 0.0467 | 0.0448 | +0.0019 [-0.0014, +0.0052] | inconclusive |
| band, frozen | logistic | 78 | 0.4401 | 0.4421 | -0.0021 [-0.0228, +0.0235] | |
| band, amended | logistic | 61 | 0.4341 | 0.4267 | +0.0074 [-0.0101, +0.0295] | |

"Inconclusive" here means no group's upper bound reaches its bar, so doubling the replications could not
change the outcome; the move stays out by default. Neither reading shows a gain, and the amended reading
shows no harm.

The amended reading was chosen after the probability results were seen, so it describes the pilot and does
not test it. How much the choice matters: the classification depends on two cutoffs, the truth's distance
from 0 and 1 (0.005, 0.01, 0.02) and the share of interior rows below which a cell counts (25, 50, 75
percent). Probability-group gain of the move, with bounds, and cells left out of 44:

| distance, share | probit left out | probit gain [bounds] | logistic left out | logistic gain [bounds] |
|---|---|---|---|---|
| 0.005, 25 | 11 | -0.0053 [-0.0117, -0.0001] | 11 | +0.0024 [-0.0035, +0.0094] |
| 0.01, 25 | 12 | -0.0022 [-0.0083, +0.0028] | 12 | -0.0003 [-0.0055, +0.0072] |
| 0.02, 25 | 13 | -0.0031 [-0.0080, +0.0019] | 13 | -0.0007 [-0.0060, +0.0070] |
| 0.005, 50 | 14 | -0.0013 [-0.0056, +0.0034] | 14 | +0.0065 [+0.0027, +0.0097] |
| 0.01, 50 | 16 | +0.0007 [-0.0034, +0.0050] | 17 | +0.0019 [-0.0014, +0.0052] |
| 0.02, 50 | 17 | +0.0010 [-0.0032, +0.0054] | 17 | +0.0019 [-0.0014, +0.0052] |
| 0.005, 75 | 18 | +0.0009 [-0.0033, +0.0053] | 17 | +0.0019 [-0.0014, +0.0052] |
| 0.01, 75 | 22 | +0.0007 [-0.0041, +0.0060] | 20 | +0.0021 [-0.0015, +0.0055] |
| 0.02, 75 | 25 | +0.0004 [-0.0036, +0.0049] | 24 | +0.0004 [-0.0035, +0.0041] |

No setting shows harm, and none shows a gain near the 0.02 bar; the largest, logistic at (0.005, 50), is
+0.0065 and is the only one whose lower bound is above zero. Under probit at (0.005, 25) the upper bound on
the worsening is 0.0117 against the 0.010 margin, so the move is there neither shown harmful nor shown
harmless. The band's gain runs from -0.0189 to -0.0048 under probit and from +0.0023 to +0.0086 under
logistic across the nine settings.

Cells with 23 to 50 percent of rows inside [0.01, 0.99] sit between the extremes and are not clearly in the
move's favor. In logistic climate, wdbc and ionosphere (BART-fitted truth, 39, 23 and 50 percent) the move
covers 0.68, 0.89 and 0.90 against the shipped sampler's 0.88, 0.99 and 0.98; on the rows inside
[0.01, 0.99] alone it covers 0.92, 0.97 and 0.94 against 0.97, 0.97 and 0.97. So the amended gain depends a
little on where the line is drawn.

## The real-outcome diagnosis, corrected

The question: what share of the long reference's draws fall inside a default-length 90 percent interval
(the reference mass)? It is 0.90 if the interval reproduces the posterior. Its floor, what two halves of the
reference give each other, is a little under that.

As first computed, over all rows, under logistic the shipped sampler (equal time) scored 0.746 against the
move's 0.873, with a floor of 0.906 to 0.908; under probit 0.859 against 0.865. That looked like a large
mixing gap. It was mostly an artifact of two things.

- **Single-precision rounding on saturated rows.** The references and stored draws were single-precision
  probabilities, which round any value within about 3e-8 of 1 to exactly 1 and tiny values to 0. Under
  separation the references' mean width was exactly zero on banknote, mushroom, tictactoe and wdbc (so width
  ratios came out infinite), and 52 percent of the move's logistic interval ends on mushroom fell where the
  rounding flips a reference draw to "outside".
- **Averaging.** The mean over rows silently dropped every row where one column was undefined. Recomputed
  over all rows the review had logistic 0.746 against 0.864 and probit 0.840 against 0.875.

Read on rows whose reference interval (5th to 95th percentile) lies inside [0.01, 0.99], with each column
averaged over the rows where it is defined (52 percent of all reference rows under logistic, 49 percent under
probit), in double precision:

| link | shipped at equal time | move at default | floor | accepted references only (shipped, move) |
|---|---|---|---|---|
| logistic | 0.872 | 0.881 | 0.899 | 0.895, 0.900 |
| probit | 0.878 | 0.876 | 0.895 | 0.898, 0.898 |

The width ratio is now finite, 0.976 to 0.984 for the two readings. So at the default length both samplers
reproduce the long-run intervals on interior rows about equally well; the move is ahead by about one point
under logistic on all references and by half a point on accepted ones, and level under probit. Most of the
earlier logistic gap came from saturated rows, not mixing. The saturated rows are where the shipped sampler
differs, and there the move is the one that reaches the long-run posterior (first computation, logistic
mushroom: reference mass 0.64 against 0.24, from the review).

Caveat: only 18 of the 52 references (26 problems a link) pass acceptance, 12 under probit and 6 under
logistic. Acceptance asks that the two halves agree, that split-Rhat be at most 1.01 for k and for the mean
held-out probability and at most 1.02 at the 99th percentile of rows, and that tail effective draws be at
least 400. Every failure includes k's split-Rhat above 1.01 (20 of 26 logistic problems, 14 of 26 probit); the
separable and large datasets fail. The diagnosis is
conditional on those 26 problems, and on separable data the long runs themselves are not settled.

## Timing

The move's cost is the ratio of sampling time per iteration, move over shipped, per cell, on one x86 core,
from three alternating timed runs of 500 iterations. 210 cell-links (105 a link) were timed, 0.61 core-hours.

| problem family | probit | logistic |
|---|---|---|
| probabilities on 22 datasets | 1.12 (1.02 to 1.27) | 1.11 (1.03 to 1.16) |
| IHDP-based effects | 1.06 | 1.04 |
| ACIC 2016-based effects | 1.06 | 1.05 |
| ACIC 2019 effects | 1.06 | 1.07 |
| Twins | 1.08 | 1.09 |
| four-arm design | 1.02 | 1.05 |
| real outcomes, 22 datasets and adult at full size | 1.13 (1.30 at adult full size) | 1.11 |
| real outcomes, three causal datasets | 1.04 | 1.06 |

Every ratio is between 1.006 and 1.299; the mean is 1.097 under probit and 1.089 under logistic. The cost is
smallest on the effect designs. Giving the shipped sampler the extra time barely changes the probability
reading: the move's gain is -0.0209 at equal iterations and -0.0200 at equal time under probit, -0.0595 and
-0.0594 under logistic; at the long reading, -0.0118 and -0.0637.

## What the pilot cannot settle

- **Whether slow k mixing harms any estimand users report.** Every probability truth is a posterior mean,
  which no interval need cover near 0 or 1, so the probability group cannot tell a sampler that reaches the
  posterior from one that does not. It would need truths drawn from the model's own prior, or a
  trusted long-run posterior that both samplers reach.
- **The effect and false-declaration groups at the sample size that matters.** 10 replications a cell: the
  bounds on the zero-coverage group are tight (about 0.004 either side), but the any-unit-declared reading
  under probit (11 cells, 110 pairs: -0.034 to +0.045) and the four-arm design (4 cells, 40 pairs: -0.094 to
  +0.125) cannot resolve a bar of 0.05. The logistic effect gain (+0.0085) is a real but small difference, with an upper bound of 0.0140
  against the 0.02 bar.
- **Regimes not tried.** Ratio estimands at rare rates, propensity weighting, rare outcomes at large n,
  stan4bart's varying intercept and a dbartsSampler loop.
- **The amended probability reading** (above) was chosen after the probability results were seen.
- **The truths themselves.** The BART truths for adult, bank, magic and spambase had not converged in k
  (split-Rhat 1.04 to 1.19); both builds face the same fixed surface, so this adds no bias between them, but
  it is a flaw in four datasets' truths.
- **The real-outcome diagnosis** rests on 18 of 52 accepted references.

## Why the study stopped, and what replaces it

Whether k needs to move for binary fits was settled earlier, in the affirmative
([binary-k-fixed-or-drawn.md](binary-k-fixed-or-drawn.md); the maintainer restated it in dec-B462). The pilot's
question was a narrower one: whether this particular move improves what users estimate. The pilot found no such gain on effects or
false declarations and found its probability reading invalid near separation, which is where the moves
differ. Finishing the main run (about 25 to 30 more core-hours, its outcome already set by the pilot) would
have asked the probability question again with the same truths. The maintainer on 2026-10-17 stopped the
study after the pilot (dec-B462), disclosed as an unplanned stop for futility, and asked first whether slow k
mixing harms the estimands users report and where, and only then what a faster move costs. The pilot used
about 30 core-hours.

The replacing study ([binary-k-mixing-need.md](../plans/binary-k-mixing-need.md)) draws truths from the BART
prior on real covariates, fits an oracle at the true k, and compares the shipped sampler, the move and k fixed
at its posterior's quantiles against long-run references that both samplers reach, across the regimes the
pilot did not try; logistic cells count as the same cell under probit (dec-B463). It runs after the merge to
main and before 1.0-0 (dec-B462). The collapsed move is neither dropped nor adopted on this evidence.
