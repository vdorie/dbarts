# binary-k-mixing-screen: what the stored references predict for the k-mixing study

Status: PREDICTION, 2026-10-29. Written before the study's first tranche, to be checked against it.

The k-mixing study ([binary-k-mixing-need.md](../plans/binary-k-mixing-need.md)) asks whether slow mixing of the leaf-prior
scale k harms what binary-fit users report. Before any of its fits run, this screen reads the binary interval
study's 52 long reference runs (26 a link, on 23 real prediction problems and three real causal datasets),
which store k beside every draw, and writes down what the study's first tranche should find. It decides
nothing; it is there so the main run can be checked against a prediction made in advance.

## What was measured

For every quantity a reference holds whose interval sits inside the interior (a probability whose 90 percent
interval lies within [0.01, 0.99], or an effect whose interval has positive width), the reference's draws
were grouped by k:

- **Variance share of log k**: how much of the quantity's posterior variance is explained by which tenth of
  k's posterior the draw came from. Near 0 means the quantity barely depends on k.
- **Mass short, one chain at a k quantile**: take only the draws whose k lies near k's 10th (or 90th)
  percentile, as a chain stuck there would, form their 90 percent interval, and measure how much of the
  whole posterior it holds. Reported as 0.90 minus that mass, so a positive number is an interval that is
  too narrow or off center.
- **Stuck-run proxy**: the study's "pooled stuck run" imitated inside the reference. Four chains, each
  stuck near one of seven quantiles of k (5th to 95th, chosen as often as the share of k's posterior each
  stands for), 500 draws each; their pooled 90 percent interval's shortfall from 0.90, averaged over 200
  such runs and over the problem's interior quantities.

The study calls a cell "no problem, bounded" when two things hold at the 0.02 bar: the pooled stuck run's
shortfall, and k's share of the start-up transient. Only the first can be read from stored references, so the
prediction below concerns that half.

## The prediction

Every regime the pilot problems stand for meets the stuck-run half of "bounded" with room: the worst proxy is
0.010 (kyphosis under logistic) against a bar of 0.02, and most are below 0.005. One chain held still at an
extreme of k does lose real coverage (up to 0.15 at k's 90th percentile for the smallest datasets), but runs
whose chains sit at different values of k, as stationary chains do, average that away.

So tranche 1 is predicted to find:

- **Interior probabilities, moderate control (S6), rare outcomes (S2, S3)**: stuck-run shortfall well
  inside the bar under both links.
- **Small samples with strong signal (S4)**: inside the bar under probit; under logistic the worst problem
  sits at half the bar, so with 12 references its upper bound may not clear 0.02. Predicted bounded under
  probit, unresolved or barely bounded under logistic.
- **Effects (S7, S8, S9)**: average and individual effects well inside the bar under both links.
- **Near separation (S1)**: no prediction. Its four datasets' references hold no interior quantity at all
  (every probability's interval reaches past 0.01 or 0.99), so the interior classes the study decides on
  were not present to read. The study reads these as tail rows, a tier 3 class.
- **Not covered by any stored reference**: many predictors (S5), stan4bart (E1, E2), the sampler loop (E3),
  the multiple-treatments design (D1) and the weighted effect (D2). No prediction for them.

Whether "bounded" is reached also depends on k's share of the transient, which the screen cannot see, and
on the references passing acceptance.

## By regime

The worst and mean stuck-run proxy over the problems standing for each regime. The prediction reads
"bounded" when the worst problem is below 0.01, half the bar, leaving room for the noise of 12 references.

| regime | link | class | problems | worst proxy | mean proxy | prediction |
|---|---|---|---|---|---|---|
| S2 | probit | interior probabilities | 3 | 0.0067 | 0.0043 | bounded |
| S2 | logistic | interior probabilities | 3 | 0.0044 | 0.0038 | bounded |
| S3 | probit | interior probabilities | 1 | 0.0021 | 0.0021 | bounded |
| S3 | logistic | interior probabilities | 1 | 0.0013 | 0.0013 | bounded |
| S4 | probit | interior probabilities | 5 | 0.0079 | 0.0060 | bounded |
| S4 | logistic | interior probabilities | 5 | 0.0102 | 0.0076 | at the margin |
| S6 | probit | interior probabilities | 2 | 0.0067 | 0.0043 | bounded |
| S6 | logistic | interior probabilities | 2 | 0.0040 | 0.0035 | bounded |
| effects (S7, S8, S9) | probit | individual effects | 3 | 0.0030 | 0.0016 | bounded |
| effects (S7, S8, S9) | probit | average effect | 3 | 0.0015 | -0.0016 | bounded |
| effects (S7, S8, S9) | logistic | individual effects | 3 | 0.0038 | 0.0027 | bounded |
| effects (S7, S8, S9) | logistic | average effect | 3 | 0.0054 | 0.0032 | bounded |

## By problem

The adult rows are the 4,000-row training subsample; adult-full is the full dataset, read beside it. Problems
missing from the table hold no interior quantity in either link (banknote, climate, ionosphere, mushroom,
sonar, tictactoe, wdbc; biopsy has one under logistic only).

| problem | link | class | interior quantities | pilot acceptance | variance share of log k (mean, max) | mass short, one chain at k q10 | mass short, one chain at k q90 | stuck-run proxy |
|---|---|---|---|---|---|---|---|---|
| actg175 | logistic | average effect | 1 | passed | 0.007, 0.007 | 0.003 | 0.011 | 0.0005 |
| actg175 | probit | average effect | 1 | passed | 0.008, 0.008 | -0.015 | 0.004 | -0.0035 |
| lalonde | logistic | average effect | 1 | failed | 0.005, 0.005 | -0.021 | 0.035 | 0.0036 |
| lalonde | probit | average effect | 1 | passed | 0.006, 0.006 | -0.013 | 0.030 | 0.0015 |
| lindner | logistic | average effect | 1 | failed | 0.028, 0.028 | 0.011 | 0.012 | 0.0054 |
| lindner | probit | average effect | 1 | failed | 0.021, 0.021 | -0.002 | 0.011 | -0.0028 |
| adult | logistic | interior probabilities | 531 | failed | 0.005, 0.119 | -0.009 | 0.019 | 0.0013 |
| adult | probit | interior probabilities | 496 | failed | 0.006, 0.105 | -0.007 | 0.017 | 0.0021 |
| adult-full | logistic | interior probabilities | 1151 | failed | 0.022, 0.128 | 0.024 | 0.011 | 0.0036 |
| adult-full | probit | interior probabilities | 1089 | failed | 0.015, 0.140 | 0.017 | 0.007 | 0.0025 |
| bank | logistic | interior probabilities | 360 | failed | 0.005, 0.031 | -0.009 | 0.018 | 0.0015 |
| bank | probit | interior probabilities | 330 | failed | 0.007, 0.046 | -0.008 | 0.021 | 0.0023 |
| biopsy | logistic | interior probabilities | 1 | failed | 0.008, 0.008 | -0.065 | 0.151 | 0.0086 |
| birthwt | logistic | interior probabilities | 38 | passed | 0.044, 0.135 | -0.030 | 0.103 | 0.0098 |
| birthwt | probit | interior probabilities | 37 | passed | 0.045, 0.146 | -0.028 | 0.109 | 0.0064 |
| cleveland | logistic | interior probabilities | 26 | failed | 0.013, 0.092 | -0.044 | 0.084 | 0.0038 |
| cleveland | probit | interior probabilities | 26 | passed | 0.013, 0.096 | -0.043 | 0.081 | 0.0018 |
| creditapproval | logistic | interior probabilities | 50 | failed | 0.017, 0.079 | -0.035 | 0.080 | 0.0044 |
| creditapproval | probit | interior probabilities | 42 | passed | 0.018, 0.112 | -0.033 | 0.070 | 0.0042 |
| german | logistic | interior probabilities | 188 | passed | 0.011, 0.079 | -0.022 | 0.046 | 0.0040 |
| german | probit | interior probabilities | 167 | passed | 0.008, 0.045 | -0.022 | 0.042 | 0.0019 |
| haberman | logistic | interior probabilities | 62 | passed | 0.031, 0.177 | -0.028 | 0.090 | 0.0069 |
| haberman | probit | interior probabilities | 62 | passed | 0.030, 0.157 | -0.029 | 0.088 | 0.0060 |
| infert | logistic | interior probabilities | 50 | passed | 0.050, 0.220 | -0.015 | 0.095 | 0.0072 |
| infert | probit | interior probabilities | 50 | passed | 0.043, 0.192 | -0.023 | 0.088 | 0.0079 |
| kyphosis | logistic | interior probabilities | 11 | passed | 0.027, 0.097 | -0.056 | 0.150 | 0.0102 |
| kyphosis | probit | interior probabilities | 9 | passed | 0.032, 0.100 | -0.059 | 0.140 | 0.0078 |
| magic | logistic | interior probabilities | 610 | failed | 0.004, 0.029 | -0.008 | 0.017 | 0.0019 |
| magic | probit | interior probabilities | 543 | failed | 0.008, 0.077 | -0.006 | 0.017 | 0.0013 |
| pbc | logistic | interior probabilities | 50 | failed | 0.013, 0.105 | -0.041 | 0.093 | 0.0044 |
| pbc | probit | interior probabilities | 48 | passed | 0.009, 0.038 | -0.045 | 0.083 | 0.0052 |
| pima | logistic | interior probabilities | 89 | failed | 0.014, 0.077 | -0.032 | 0.071 | 0.0031 |
| pima | probit | interior probabilities | 77 | passed | 0.013, 0.081 | -0.032 | 0.066 | 0.0067 |
| spambase | logistic | interior probabilities | 200 | failed | 0.008, 0.054 | -0.012 | 0.024 | 0.0019 |
| spambase | probit | interior probabilities | 181 | failed | 0.011, 0.078 | -0.005 | 0.021 | 0.0016 |
| transfusion | logistic | interior probabilities | 148 | failed | 0.016, 0.085 | -0.021 | 0.048 | 0.0015 |
| transfusion | probit | interior probabilities | 143 | passed | 0.015, 0.080 | -0.020 | 0.053 | 0.0039 |
| actg175 | logistic | individual effects | 1054 | passed | 0.006, 0.036 | -0.021 | 0.038 | 0.0013 |
| actg175 | probit | individual effects | 1054 | passed | 0.007, 0.033 | -0.027 | 0.036 | -0.0004 |
| lalonde | logistic | individual effects | 614 | failed | 0.012, 0.074 | -0.033 | 0.064 | 0.0038 |
| lalonde | probit | individual effects | 614 | passed | 0.014, 0.079 | -0.036 | 0.067 | 0.0021 |
| lindner | logistic | individual effects | 996 | failed | 0.012, 0.128 | -0.026 | 0.041 | 0.0029 |
| lindner | probit | individual effects | 996 | failed | 0.012, 0.131 | -0.028 | 0.049 | 0.0030 |

## Caveats

- Within a reference, keeping only the draws near one value of k approximates a chain run with k fixed
  there, but it is not one: those draws come from chains that moved through k, so their trees are mixed as
  the reference's are. A fixed-k fit of default length also carries its own transient.
- 34 of the 52 references failed the pilot's acceptance (its halves disagreed or did not converge), most of
  them on near-separated or large problems; their draws are read as they are.
- These are real outcomes with one fit each, not draws from the model's prior, and the regimes they stand
  for are matched by dataset, not by the study's selection on the data.
- The proxies are means over interior quantities; the study reads the same means, so a few quantities with
  large shortfalls can hide inside a small mean in both.
