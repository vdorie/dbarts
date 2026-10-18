# binary-interval-study: does the collapsed scale move make binary intervals more accurate at equal running time?

Status: PLANNED 2026-10-17, revised from its blind critique (scratch/bisc-crit.md) and the orchestrator's
calls on it. The bar is ruled (dec-B455): calibration, distance from nominal on either side; real outcomes
against a reference are a diagnosis outside the rule (dec-B456). Nothing run
beyond the planning probes listed at the end. The maintainer placed their multiple-comparisons problem in a
multiple-treatments scenario; the design A6 answers it ([Multiple comparisons](#multiple-comparisons)). The
run waits on nothing from the maintainer.

agent: one opus implementer writes the harness, runs it on the bench box and writes the finding.
rng: NEUTRAL. No package code changes; a benchmark harness and two private builds of the move's branch.
window: before 1.0-0. It decides both halves of the collapsed move, probit and logistic, which the
maintainer holds (2026-10-17).
budget: about 49 core-hours expected; the decision path, doubled replications included, fits in 74; the
stop is at 75 ([Cost](#cost)). Harness about 1,500 lines of R.

## Summary for the maintainer

You asked for proof that the collapsed scale move helps with what users estimate, not with k: intervals
that are too narrow, and a multiple-comparisons problem, on real binary problems, at equal running time.
What exists is one number: logistic fits on simulated 500-row near-separable data covered the true
probability 79 percent of the time without the move and 99 percent with it. That was at 64,000 draws, not
the default, and on a probit truth that a logistic fit's tails cannot follow. Under your ruling that
calibration is the goal, it is a move from 11 points under nominal to 9 points over: barely a gain.

The study compares the shipped sampler with the same sampler plus the move, for probit and for logistic.

- **Real covariates, made-up outcomes, so the truth is known.**
  - The 22 real binary datasets of the k study, each outcome replaced by draws from a surface fitted to the
    real outcome. Two such surfaces per dataset: one from BART, one from a smooth additive model, so the
    answer does not rest on BART's own prior.
  - IHDP's infant-development covariates with treatment effects built on them, with a binary outcome.
  - The ACIC 2016 competition's 58 real covariates and generated surfaces turned into binary outcomes,
    including a version with 250 covariates, effects confined to a subgroup, an outcome near 1 percent and
    rare outcomes where treated and control units barely overlap.
  - The ACIC 2019 competition's eight binary practice datasets, which come with their true probabilities.
  - Twins: real infant mortality where both twins' outcomes are observed, so the effect of being the heavier
    twin is known on real outcomes, for averages and subgroups.
  - Designs where every true effect is zero, and one where effects sit in a subgroup only, which is where
    a multiple-comparisons problem shows.
  - The case you named: several treatments, here four arms on the ACIC 2016 covariates, some pairs of arms
    with real effects and some with none, one arm with weak overlap, and a cell where no pair differs.
- **Real outcomes, truth unknown.** The same 22 datasets and the adult census data at full size with their
  real outcomes, and three real causal datasets (a heart-procedure study with rare deaths, an HIV trial,
  the LaLonde job-training data). Each gets one very long run, half from each sampler, as the gold
  standard; the question is whether default-length runs reproduce its intervals.
- Each fit is read at the default length, at the length the shipped sampler can run in the time the move
  takes (the move costs 14 to 28 percent more per iteration at the default 75 trees), and at a longer
  length.

What is measured: whether 90 and 95 percent intervals contain the truth as often as they claim, neither more
nor less often, for the predicted probability and each unit's treatment effect; their width; whether all
units' intervals hold together; how often a unit with no effect is declared to have one, and how often a
declared effect has the wrong sign; and the error of the posterior mean. The study also asks what goes wrong
that the move cannot touch: whether the exact posterior's own intervals miss (the model or prior), whether
the trees mix too slowly even with the move, and where (rare outcomes, near 0 or 1, no comparable treated
units, many covariates).

The rule, fixed before any result: a link's move is kept if, at equal running time, it brings coverage at
least 2 points closer to nominal on one of three groups (probabilities with known truth, effects with known
truth, false declarations of effects) and is shown to do no harm
on any. It is left out if harm is shown. Anything else is inconclusive: the replications are doubled once,
and a half still inconclusive stays out by default, the numbers going to you.

Cost: about 49 core-hours, about 14 hours of the box on its four cores, less if the halfway reading settles
it; up to 74 if every replication has to be doubled.

The multiple-comparisons problem you placed in a multiple-treatments scenario is a design of its own, read for
false declarations among the pairs of arms with no effect, for all contrasts within a unit at once and across
units, and for wrong signs. It is part of the decision, and harm on it alone leaves the move out. The exact
case your colleagues saw, if they report it, can be scored on the stored draws or run as a new cell after the
study.

## What is known

- Binary fits draw k under chi(1.5, 2) by default (dec-B420). Probit carries a cheap rescaling step that
  costs 0 to 3 percent of an iteration (dec-B371, dec-B397); at 500 rows and 75 trees it did not converge
  k ([binary-k-fixed-or-drawn.md](../design/binary-k-fixed-or-drawn.md): split-Rhat 1.12 to 1.21 at
  64,000 draws). Logistic has no cheap step, its latents being Polya-Gamma precisions
  ([k-mixing-pg-families.md](../design/k-mixing-pg-families.md)).
- The collapsed move is built on the branch wt/collapsed-move-plan, tip 207f740e (based on bartcore
  3ff28d4c), not on bartcore: the leaf and k helper, the logistic half, the probit half, its tinytest and
  the exact gate's logistic arms. Its early measurements (scratch reports csm, csmr, csmx, csml; the
  count-weight question settled in csmq) give:
  - per-iteration cost on the x86 box, with the move over without: probit 1.17 at 500 rows and 75 trees,
    1.25 at 2,000, 1.28 at 5,000, 1.08 to 1.16 at 200 trees; logistic 1.14, 1.17, 1.15 at 75 trees,
    1.07 to 1.11 at 200;
  - at the k study's near-separable 500-row cells (8 chains of 8,000 after 8,000): logistic's k converged
    in every fit with the move (none without); probit's in 73 percent (4 percent without); the held-out
    probability's split-Rhat above 1.05 in 40 percent of probit fits with the move, 65 without;
  - the one interval number, from the same runs: logistic's median 90 percent coverage of the true
    probability 0.79 without the move, 0.99 with it, a probit truth throughout;
  - on large fits (10,000 and 50,000 rows) the gains in k's autocorrelation carry a relative error near
    one half.
- The builds those reports used: `csm`, a plain install of the branch tip; `csm-off`, the same source with
  the engine's test toggle for the move, a default of true, flipped off by a one-line sed (the patch file
  kept beside it does not apply); on the box, two exports of 207f740e built the same way. csml found the
  flipped build bit for bit equal to bartcore d5b637ee on arm64, on seeded probit and logistic fits with an
  offset, counts and a mask; the source between d5b637ee and 3ff28d4c changes the variance forest only. The
  critique's probe found it bit for bit equal to bartcore 92ae2714 on arm64 too, under both links, with
  numeric predictors and with a 4-level factor, test rows and 2 chains; the move build differed on all four.
  So both arms are built from the branch, with no rebase (the orchestrator's call), and stage 1 repeats the
  check on x86.
- The k study's real-data harness, datasets and loaders are benchmarks/R/binary-hyperprior.R and
  benchmarks/R/uci-binary.R: equal-tailed intervals of the per-draw probability, coverage distance against
  0.90, interior rows split out. [benchmark-surfaces.md](../design/benchmark-surfaces.md) reads coverage
  against a well-mixed reference arm where the shipped configuration over-covers; this study reports that
  reading wherever a reference exists, beside the rule's reading against nominal.

## Problems

### A. Truth known: real covariates, generated outcomes

**A1. The 22 real datasets with fitted surfaces as truth.** Two truths per dataset and link, each fitted to
the real outcome on the whole dataset (a fixed 5,000-row subsample for the five large ones):
- **BART truth**: the posterior mean probability of the move build at reference length (16 chains of
  5,000 burn-in and 15,000 iterations). The move build and the long run keep a short chain's k bias out of
  the truth, which would favor the shipped arm.
- **Additive truth**: a penalized additive model of that link (mgcv): smooth main effects for numeric
  predictors with ten or more distinct values, the rest as terms of their own, and tensor-product
  interactions among the three numeric predictors with the largest main-effect deviance reduction; on
  datasets of more than 20 predictors, smooths on the 20 with the largest univariate deviance reduction and
  the rest linear; a ridge penalty on the parametric terms where the data separate; probabilities held in
  [1e-4, 1 - 1e-4]. Out of BART's class, and not shaped by its prior.

Each replication draws a new 80/20 split (4,000 and 1,000 rows for the large five) and a new outcome
y ~ Bernoulli(p(x)) on every row, and scores the held-out intervals against p. Each link is scored on its
own link's truths; a cross-link reading (each link on the other link's BART truth, 3 replications, shipped
build only) is diagnosis. 44 cells a link, 16 replications each.

**A2. IHDP with a binary outcome.** The 747 children of the IHDP construction after its deletion rule
(every treated child of a nonwhite mother dropped), 25 standardized covariates, the real treatment. Latent
surfaces follow npci's settings A and B (npci is the maintainer's implementation; Hill 2011 labels them
differently, its exponential surface being B): A, a nonlinear control surface exp((x + 0.5) b) against a
linear treated surface x b - w, b sparse; B, main effects, quadratics and interactions with sparse
coefficients for both arms. Each surface is standardized on the controls and mapped to p_z(x) = F(a eta_z(x)
+ c), F the fitted link, c set for a control base rate of 0.3, a = 1 (moderate) or 2.5 (near-separated).
The individual effect is p_1(x) - p_0(x); null cells set p_1 = p_0. Units whose mother is nonwhite have no
treated comparison and are read separately. Cells: A moderate, A near-separated, B moderate, null moderate,
null near-separated; 30 replications each.

**A3. ACIC 2016 covariates with a binary outcome.** The competition's 4,802 rows and 58 real covariates
(Collaborative Perinatal Project: 3 categorical, 5 binary, 27 counts, 23 continuous), its treatment
assignment and its generated surfaces, mapped to probabilities as in A2. Settings: row 1 (linear surface,
high heterogeneity, weak overlap: propensities down to 1e-10), row 31 (exponential surface, medium
heterogeneity, full overlap, low alignment between what drives treatment and what drives the outcome), row
30 (step surface, high heterogeneity, full overlap), and a null built on row 31's control surface. Cells:
- at 1,000 units, a new subsample each replication, 30 replications: the three settings and the null at
  both signals, base rate 0.3 (8 cells); row 31 and its null at a base rate of 0.05 (2); row 1 at a base
  rate of 0.05, a rare outcome under weak overlap (1); sparse heterogeneity, row 31's control surface with
  an effect of one latent standard deviation in a subgroup of about 15 percent of units defined by two
  covariates and none elsewhere (1);
- high-dimensional, 1,000 units, 30 replications: row 31 and its null with the 58 covariates and 192
  permuted copies of them (each a real column with its rows permuted independently), p = 250 (2);
- at all 4,802 units, 16 replications: row 31 and its null, near-separated (2); row 31 at a base rate of
  0.01, about 48 events, moderate (1).

17 cells.

**A4. ACIC 2019's binary practice datasets.** Four low-dimensional (500 to 625 units, 15 to 22 covariates)
and four high-dimensional (2,000 units, 185 covariates), each with its true probability under both arms for
every unit; two of the low-dimensional ones have control probabilities from 0 to 0.96 and 0.999. A
replication draws a new outcome from the given probabilities, covariates and treatment fixed. 8 cells, 20
replications each.

**A5. Twins.** The 11,984 same-sex twin pairs both under 2,000 g of the CEVAE data (Louizos et al. 2017),
one-year mortality 0.19 for the lighter twin and 0.16 for the heavier, 10 percent of pairs discordant, 46
covariates with missing values. The heavier twin is treated; each replication subsamples 2,000 pairs and
assigns treatment by Louizos et al.'s covariate-driven rule, observing one twin. Both outcomes are real, so
the sample average effect and the averages in the ten gestational-age groups are known exactly; a unit's
effect is a realized 0 or 1 difference, not a probability, so unit intervals are not scored. The
counterfactual twin scores the predicted counterfactual probability (log score, calibration). 1 cell, 20
replications.

**A6. Multiple treatments, the maintainer's case.** The ACIC 2016 covariates, 2,000 units subsampled each
replication, four arms (0 to 3) assigned by a multinomial logistic model on the covariates, arm 3's
propensity driven toward zero over part of the covariate space so that it has weak overlap. Latent surfaces,
mapped to probabilities as in A2 at a base rate of 0.3: arm 0 row 31's control surface; arm 1 that surface plus
row 31's heterogeneous effect; arm 2 the same as arm 1 plus a second heterogeneous effect on other
covariates; arm 3 the same as arm 0. Of the six pairwise contrasts, 1 against 0, 2 against 0, 2 against 1, 3
against 1 and 3 against 2 are heterogeneous and 3 against 0 is exactly zero for every unit; a second
arrangement makes 2 against 1 zero instead (arm 2 equal to arm 1), so two of six pairs are null; the null
cell gives every arm arm 0's surface, so all six are. Cells: the first arrangement moderate and
near-separated, the second moderate, the null moderate; 30 replications each. The fit is bart() with the arm
as a factor predictor beside the covariates and the estimated generalized propensity scores, predicting
every unit under all four arms. A forest multiplied by the arm, the multi-forest route, does not exist here:
push 2 of forest-kind-by-class refuses a forest multiplied by a factor of three or more levels in 1.0
(dec-B281, dec-B282). 4 cells.

### B. Truth unknown: real outcomes against a gold-standard run

**B1. The 22 real datasets with their real outcomes**, one fixed split each (the k study's first), and
**adult at full size** (39,073 training rows, 2,000 fixed held-out rows scored), 12 seeds a dataset (adult
8). The
seeds differ only in the chain, so their spread is Monte Carlo error and the gap to the reference is what
poor mixing costs.

**B2. Three real causal datasets with binary outcomes**, individual effects on every unit, 15 seeds each:
- lindner (twang): 996 coronary-intervention patients, abciximab or not, six-month survival, 26 deaths.
- ACTG 175 (BART): zidovudine against zidovudine with didanosine, 1,054 patients, the event indicator (about
  27 percent), baseline covariates only. Randomized.
- lalonde (twang): 614 men, job training, earnings in 1978 above zero (77 percent). Weak overlap.

### Data and how it is obtained

| source | what | how | status |
|---|---|---|---|
| six R datasets | pima, biopsy, infert, kyphosis, pbc, birthwt | MASS, rpart, survival, datasets | offline |
| sixteen UCI datasets | adult to sonar | benchmarks/R/uci-binary.R, sha256-checked | cached on the laptop; copied to the box |
| IHDP | 747 units after deletion | ihdp.RData in the maintainer's npci checkout | offline, copied |
| ACIC 2016 | covariates, 77 settings, generator | aciccomp2016 from the maintainer's aciccomp checkout (GitHub, not CRAN) | offline; installed in planning |
| lindner, lalonde; ACTG 175 | real causal | twang; BART | offline |
| Twins | three files of 71,345 pairs | the CEVAE repository | fetched in planning |
| ACIC 2019 | eight binary practice datasets with truth (of 16; the rest continuous) | the challenge's archive of practice files | fetched in planning |

Every download goes into its own new, empty directory, is checked against a recorded sha256, and is read by
R running from outside that directory; scripts live elsewhere. In planning: scratch/bis-data/twins (four
files) and scratch/bis-data/acic2019 (the low- and high-dimensional practice archives, each extracted into
a directory of its own). The challenge's 3,200 scored datasets carry no unit truth and were not fetched.

Duplicate covariate patterns matter for a drawn k: rows that repeat one outcome pull the posterior toward
interpolation (csmq). The share of rows in repeated patterns with one outcome is 0.15 in transfusion, 0.12
in spambase, 0.11 in haberman, under 0.03 elsewhere among the UCI sixteen; results are split by it.

## Arms

Two builds of the branch tip, each in its own private library:
- **shipped**: the toggle flipped off. Probit runs its rescaling step, logistic nothing: the sampler bartcore
  ships.
- **move**: the tip as is. Probit runs the step and the collapsed move, logistic the collapsed move.

Package defaults otherwise: 75 trees, chi(1.5, 2), 500 burn-in, and `n.threads = 1` set explicitly on every
fit (bart() otherwise takes as many threads as chains). Fits of individual effects use bartCause's 10 chains
and its default response model, which carries an estimated propensity score: one propensity per replication
(per cell where the covariates and treatment are fixed), the posterior mean of a 10-chain BART probit fit of
the treatment on the covariates by the shipped build, passed to both arms as a covariate, so the move
touches the outcome fit only (the orchestrator's call). A6 has no bartCause fit (bartCause takes one binary
treatment); it uses the same 10 chains, and its generalized propensity is the posterior mean of a 4-chain
multinomial fit by the shipped build, one per replication, shared by both arms and both links. The other fits
use bart()'s 4 chains.

Each fit runs past the default and is read at several lengths: a longer run's first draws are bit for bit a
shorter run's with the same seed (checked in planning under both links). A 4-chain fit runs 500 burn-in and
2,000 draws a chain, a 10-chain fit 500 and 1,000. Readings:

- **default**: the first 500 draws after the 500 burn-in.
- **equal time**: the shipped sampler given the move's running time. r is the ratio of sampling time, the
  move build's over the shipped build's, per cell, from the timing pass: a sampler built by bart() with
  `samplerOnly = TRUE`, warmed, its run of 500 iterations timed alone, three alternating repetitions per
  build, four jobs at once as in the main run. Setup, prediction after the run and R's handling of the draws
  are left out, which would otherwise short-change the shipped arm. (A sampler-only run does not draw the
  same chain as bart() with the same seed, checked in planning, so it times and never fits.) The shipped
  sampler then has 1,000 r iterations a chain, spent either as more draws (500 burn-in, 1,000 r - 500 draws)
  or as more burn-in (1,000 r - 500 burn-in, 500 draws), both read from the same run. Which one is chosen on
  the pilot alone, by pooled coverage distance, and frozen for the rest (the orchestrator's call). For the
  simultaneous measures the more-draws reading is thinned to 500 draws a chain, the move's count, since a
  band built from order statistics widens with the number of draws.
- **long**: all draws of the run, for both builds: what running longer buys.

The equal-sweeps comparison is reported and does not decide.

## Reference runs

A reference is 16 chains, 8 of each build, of 5,000 burn-in and 15,000 iterations at thin 15 (16,000 kept
draws). Both kernels target the same posterior, so a reference made of one kernel would favor it wherever it
had not fully mixed; made of both, the halves check each other (the orchestrator's call).

Accepted when the two halves agree on every scored quantity (each row's, or under B2 each unit's effect's,
posterior mean and 5th and 95th percentiles within 4.5 combined Monte Carlo standard errors), rank-normalized
split-Rhat over all 16 chains is at most 1.01 for k and the mean held-out probability and at most 1.02 at the
99th percentile of rows or units, and tail effective draws at the 5th and 95th percentiles are at least 400.
A cross-check that costs nothing: the pooled long readings of the seeds (8 to 15 seeds a build, 4 or 10 chains
each, started apart) must agree with the reference by the same test.

A problem whose reference fails stays in, scored against the pooled 16 chains and flagged; R is reported with
and without the flagged problems. Agreement within one kernel's chains cannot
see stickiness that every chain shares; a failing problem is itself a finding about tree-level mixing.

References run on every B problem (the diagnosis R reads them) and, for the diagnosis, on A cells that qualify
([Stops](#stops)).

## Measures

All intervals equal-tailed, of the per-draw quantity, at 90 and 95 percent, transformed per draw and averaged
last.

Truth known (A):
- **Pointwise coverage and width** of the probability f(x) at held-out rows (A1) and of each unit's effect
  p_1 - p_0 (A2 to A4); of each of A6's six pairwise contrasts per unit and of each average contrast, with
  their calibration (coverage at 50, 80, 90 and 95 percent); of the ten gestational-age group averages (A5).
  Split by true probability inside
  [0.01, 0.99] or not, and in A2 by units with and without treated comparison.
- **Interval score** at 90 percent (Gneiting and Raftery 2007): width plus 20 times the distance by which the
  truth falls outside. It rewards narrow intervals only when they still cover, so it guards width.
- **Simultaneous coverage**: one band over all held-out rows or units (in A6, one band over a unit's six
  contrasts, and one over all units' contrasts), the rank-based band of Besag et al.
  (1995) at a common tail level, built on the latent scale for probabilities; the share of replications
  whose whole truth vector lies inside. Reported with the band's widening over the pointwise intervals and
  how often it collapses to the range of the draws.
- **RMSE** of the posterior mean against the truth on the probability scale.
- **False declarations**: a unit is declared when its 95 percent effect interval excludes zero. At units whose
  true effect is zero (the null cells and the sparse cell's 85 percent), the coverage of zero, nominal 0.95;
  the share of replications with any unit declared, pointwise and by the band. Among declared units elsewhere,
  the share with the wrong sign.
- **Average effect**: coverage of the sample average effect's interval.

Truth unknown (B):
- **Reference mass**: the share of the reference's draws inside the arm's interval, per row or unit, 0.90 for a
  90 percent interval that reproduces the posterior; its Monte Carlo floor from splitting the reference in two.
- **Width ratio** to the reference's interval, the interval score averaged over the reference's draws, and the
  posterior mean's error in reference standard deviations.
- **Simultaneous reference mass**: the share of reference draws entirely inside the arm's band, with its
  floor.

Mixing, every fit: split-Rhat and effective draws of k and of the mean held-out probability; sampling time.

## What goes wrong beyond the move

The maintainer reports trouble with binary outcomes on real data. These readings go to the report and TODO
whatever the rule decides:

- **The posterior itself misses**: the reference's own coverage of the truth, pooled by cell type, at least 3
  points from nominal with its bound excluding zero. That is the model or the prior; no sampler fixes it.
  Types: multiple treatments, A1 near-deterministic (most truth outside [0.01, 0.99]), A1 other, effects moderate, effects
  near-separated, rare outcomes, weak overlap or extrapolation, high-dimensional, null. Pooled, because one
  cell's few references cannot resolve 3 points.
- **Mixing the move does not reach**: the move's default reading at least 3 points of reference mass from the
  reference while k has converged: the trees.
- **Mixing the move does reach**: the shipped reading off and the move's not.
- **The link's tails**: the cross-link reading of A1.
- **Rare outcomes, extrapolation, many covariates, duplicated rows**: the A3 rare, high-dimensional and
  weak-overlap cells, lindner, IHDP's units without treated comparison, lalonde, the duplicate split.

## Decision rule

Fixed before any result. Per link, at equal running time: the move's default reading against the shipped
sampler's equal-time reading.

**The gap measure (dec-B455).** A cell's gap is its distance from nominal coverage, on either side: the
maintainer ruled for calibration ("I don't want to over-optimize and just have 100% coverage."). Taking the
absolute value of a noisy coverage inflates every cell's distance by about one cell standard error, more near
nominal, which would penalize the arm that is closer. Each cell's distance is therefore corrected for that
noise: with s the cell's standard error from its replications, the cell's distance is the d >= 0 whose folded
normal mean, s sqrt(2/pi) exp(-d^2 / (2 s^2)) + d (1 - 2 Phi(-d / s)), equals the observed |coverage -
nominal|, and zero when the observed distance is below s sqrt(2/pi), the mean under no true distance. The
correction vanishes once the distance is a few standard errors. Pooling a group's coverage before taking the
distance was not chosen: it lets cells over and under nominal cancel. A group's gap is the mean of its cells'
corrected distances.

Three gain groups, the gap of 90 percent intervals (95 percent for false declarations):
- **F**: held-out probabilities, truth known (A1's 44 cells);
- **T**: effects, truth known (A2 to A4's 15 non-null cells, A6's non-null contrasts in its three non-null
  cells, and A5's group averages);
- **M**: multiple comparisons, three measures, each with its own bar ([Multiple comparisons](#multiple-comparisons)):
  - M1, the coverage of zero wherever the true effect is zero (the 7 null cells, the sparse cell's 85
    percent, A6's null contrasts), its gap from 0.95, bar 0.02;
  - M2, the family-wise rate, the share of a cell's replications in which the 95 percent band declares any
    null unit or contrast (the null cells, A6's three cells with null pairs), its gap from 0.05, bar 0.05;
  - M3, the share of declared effects or contrasts whose true value has the other sign, a reduction, bar
    0.01.

Guards, each with a margin:

| guard | margin | planned standard error |
|---|---|---|
| F, T, M1 worse | 0.010, 0.015, 0.015 | 0.0019, 0.0038, 0.0032 |
| M2 worse, M3 worse | 0.075, 0.010 | 0.018, from the pilot |
| A6 alone, the maintainer's case: M1, M2 worse | 0.025, 0.11 | 0.0055, 0.029 |
| simultaneous coverage gap, pooled over the A cells | 0.035 | 0.0097 |
| average-effect coverage gap | 0.040 | 0.011 |
| pooled RMSE ratio, pooled interval-score ratio, over the A cells | 1.02 | |
| pairs skipped, the move's arm over the shipped arm | 0.5 percent of pairs | |

Each margin is at least 3.3 planned standard errors (M3's spread is unknown until the pilot); after the pilot each becomes the larger of the planned
margin and 3.3 times the pilot's standard error, then frozen.

Differences are paired by replication or seed, their bounds from a paired bootstrap of 2,000 resamples over
replications within cells. R, a diagnosis (dec-B456), is read by the same bootstrap over seeds and,
independently, over the reference's chains within each build, since the seeds share one reference and one
split; it is conditional on its 26 datasets.

Three looks, interim, final and doubled, spend one-sided alpha 0.005, 0.02 and 0.025, so every bound below is
at the look's level.

- **Include the move** for a link when a gain measure improves by at least its bar (0.02, or M2's and M3's
  own), its lower bound above zero, and
  every guard is shown harmless: the upper bound on its worsening below its margin. The gain groups are F, T
  and M (dec-B456).
- **Leave it out** when harm is shown on any guard: the lower bound on its worsening past its margin.
- **Otherwise inconclusive.** The undecided groups' and guards' cells are doubled once (budgeted), unless no
  gain measure's upper bound reaches its bar (0.02 for F, T and M1, 0.05 for M2, 0.01 for M3), in which
  case doubling cannot change the outcome and is skipped. Only A cells are doubled; B is diagnosis.
  Still inconclusive, the half is left out by default, both halves if both are inconclusive: neither has
  shown the interval gain asked for. The numbers go to the maintainer with that recommendation, the
  k-mixing evidence (csml) and the cost beside them.

95 percent intervals are read beside the 90 percent ones (and 90 beside 95 for M): a group whose gain changes
sign between the two levels counts as inconclusive for gain.

Outcomes:
- Probit included: the probit half stays, the collapsed-move plan's large-n options still to weigh. Left out:
  it is dropped and probit keeps its rescaling step alone.
- Logistic included: the logistic half ships on this study's evidence, a new ruling. Left out: it is dropped;
  logistic keeps a drawn k under chi(1.5, 2) with no scale move, and the report states what that costs: the
  reported k, and intervals near 0 or 1, depend on chain length.
- The diagnosis goes to the maintainer and TODO either way.

### What the counts can resolve

Planned pairs: F 704, T 692, M1 346, M2 316; A6 alone 120. With a per-replication spread of the paired coverage difference of
0.05 for F (the k study's paired standard errors imply it), 0.10 for T (one replication's units share a
surface) and 0.06 for M1, and a discordance of 0.1 between arms for M2's yes-or-no outcome, the standard errors are
those of the table. M2 can show only a large change: its bar is 0.05, against 0.02 for the rest. Inclusion needs a point estimate of at
least 0.02, so at a true gain of exactly 0.02 it succeeds about half the time; 95 percent power needs a true
gain of about 0.023 (F), 0.026 (T) and 0.026 (M1). On A6 alone a 2-point gain is under 4
standard errors, so its own reading guards and is reported, and gains count pooled into M. Correlation within a replication is handled by
resampling whole replications. The pilot (10 replications of every cell) rereads the spreads.

## Stops

- **Stage 1**, on the x86 box, before any study fit: the shipped build differs bit for bit from a build of
  bartcore 92ae2714 on seeded fits of both links along the study's own paths (numeric and factor predictors,
  missing predictor values, test rows, 4 and 10 chains, `n.threads = 1`); or the move build does not differ;
  or the prefix property fails on either build. Any of these stops the run.
- **Halfway**: with half the replications of every A cell in, the rule is read at the
  interim alpha. A decided link stops.
- **Diagnosis references** run last, one replication per A cell where either build's default reading is more
  than 3 points from nominal, plus four control cells a link chosen now (A1 pima and adult with the BART
  truth, A2 setting A moderate, A3 row 1 moderate), and a second replication per cell in a type with fewer
  than 8 cells qualifying. They are what the budget stop cuts first.
- **Skipped pairs**: a job past three times its cell's median time, or a fit with a non-finite draw, is
  skipped as a pair (both builds), tallied by the arm that failed; the tally is a guard. Skips above 1 percent
  of pairs stop the run.
- **Budget**: at 75 core-hours, stop and report what is in.

## Cost

Per-iteration times on one x86 core at 75 trees, shipped build, in milliseconds, from the k study's long
fits (probit 0.39 at 500 rows, 0.58 at 2,000, 1.19 at 4,000; logistic 0.50, 1.04, 2.06) and the move's
large-fit timings (csmx), with prediction at the counterfactual rows for effect designs:

| size class | probit | logistic | move's extra |
|---|---|---|---|
| A1, B1 small (about 450 training rows) | 0.38 | 0.50 | 17, 14 percent |
| A1, B1 large (4,000 and 1,000) | 1.19 | 2.06 | 25, 17 |
| A1 truths on whole data, small, large (5,000) | 0.45, 1.45 | 0.60, 2.50 | 17, 14; 25, 17 |
| IHDP and ACIC 2019 low-dimensional, with counterfactuals | 0.50 | 0.70 | 17, 14 |
| ACIC 2016 at 1,000 units; B2 | 0.60 | 1.10 | 17, 14 |
| high-dimensional, p 250 | 0.70 | 1.25 | 17, 14 |
| ACIC 2019 high-dimensional, 2,000 units | 1.00 | 1.90 | 25, 17 |
| Twins, 2,000 pairs | 0.90 | 1.70 | 25, 17 |
| ACIC 2016 at 4,802 units | 1.60 | 2.80 | 25, 17 |
| adult at full size | 4.50 | 9.50 | 25, 17 |

Iterations: a 4-chain fit 4 x 2,500 = 10,000; a 10-chain fit 10 x 1,500 = 15,000; a propensity fit 10 x
1,000 = 10,000; a reference or a BART truth 16 x 20,000 = 320,000. A pair costs t x iterations x (2 + extra);
a reference t x 320,000 x (1 + extra / 2); a truth t x 320,000 x (1 + extra). A planning probe (747 IHDP units
with counterfactuals, 4 chains of 1,000 iterations, the laptop under load) took 1.7 s under probit and 2.1 s
under logistic, in line.

| stage | what | probit | logistic |
|---|---|---|---|
| 1 | builds, checks, timing pass | 0.75 | 0.75 |
| 2 | A1 BART truths, 22 at reference length | 1.60 | 2.33 |
| 2 | A1 fits on BART truths, 352 pairs | 1.22 | 1.80 |
| 2 | A1 fits on additive truths, 352 pairs (the additive fits 0.03) | 1.25 | 1.83 |
| 2 | A1 cross-link, 66 shipped fits | 0.10 | 0.16 |
| 2 | A2, 150 pairs and propensities | 0.89 | 1.14 |
| 2 | A3 at 1,000 units, 360 pairs | 2.55 | 4.13 |
| 2 | A3 high-dimensional, 60 pairs | 0.50 | 0.79 |
| 2 | A3 at 4,802 units, 48 pairs | 0.93 | 1.43 |
| 2 | A4, 160 pairs | 1.44 | 2.21 |
| 2 | A5, 20 pairs | 0.22 | 0.36 |
| 2 | A6, 120 pairs (2,000 units predicted under 4 arms) and the multinomial propensities | 2.04 | 3.28 |
| 3 | B1 fits, 272 pairs (adult 8 of them) | 1.14 | 1.81 |
| 3 | B2 fits, 45 pairs | 0.25 | 0.45 |
| 3 | B references, 26 | 1.84 | 3.03 |
| | without diagnosis | 16.7 | 25.5 |
| 4 | diagnosis references, at most one per A cell (79) | up to 5.3 | up to 8.5 |
| | expected, half the diagnosis references | 19.4 | 29.7 |

Core-hours, stage 1 split evenly and A6's shared propensities (0.48) split between the links. Expected about
49 (42.2 + 13.8 / 2), about 14 hours of the box at one job a core (four jobs share memory bandwidth: 15 percent
added). With the halfway stop firing on both links, about 34. Doubling every fit on both links adds 31.2
(probit 12.2, logistic 19.0), so the decision path at its worst is 42.2 + 31.2 = 73.4; with every diagnosis
reference as well, 87. The stop at 75 covers the decision path; diagnosis references beyond it are cut. The
study runs under the bench-box grant (dec-B404), reported to the maintainer.

Memory: the largest jobs (A3 at 4,802 units at both treatment values, A6 at 2,000 units under four arms, 10
chains of 1,000 draws) hold 0.6 to 0.8 GB of draws and twice that at peak; they run at most two at a time. adult's reference scores 2,000
rows (about 0.3 GB). Disk: effect draws as single precision, thinned to 500, at the default and equal-time
readings of every effect fit, about 17 GB, so the colleagues' exact case, if reported, can be scored without
new fits.

## Scripts to write

- **benchmarks/R/binary-interval-study.R**, tracked, about 1,500 lines: the problem registry (cells and seeds),
  the generators (A1's two truths, npci's surfaces on IHDP, ACIC 2016 through aciccomp2016, the permuted
  copies, the sparse, rare and null variants, ACIC 2019's redraws, Twins' assignment and subsample), the
  propensity fit, the fit runner with its readings, the timing pass, every measure, the reference runner and
  its acceptance, and actions `jobs`, `truth`, `fit`, `reference`, `timing`, `pilot`, `interim` and
  `summary`, the last two computing the rule exactly as written here. It sources benchmarks/R/uci-binary.R;
  the six R-dataset loaders move from binary-hyperprior.R into a file both scripts source, order and seeds
  unchanged. Downloads go through a fetcher on uci-binary.R's pattern (recorded sha256, a new directory per
  source); local files (IHDP) come by environment variable.
- The box driver, in scratch: a queue of four workers, each job a pair (both builds, alternating order), jobs
  in random cell order so the halfway reading covers every cell, the largest jobs two at a time.
- The flipped build: a correct one-line patch setting the move's toggle default to false, applied to an export
  of the tip.

Outputs stay in scratch/bis; the finding goes to docs/design/binary-interval-study.md.

## Order of work

1. Write the harness; smoke it on the laptop with one replication and 50 iterations a cell, both links and
   builds.
2. On the box: both builds from 207f740e with `--preclean`, aciccomp2016, the data. Stage 1's checks, then the
   timing pass.
3. The A1 truths. The pilot: 10 replications of every cell; the equal-time choice and the margins frozen.
4. Stages 2 and 3 to the halfway point, the B references among them; the halfway reading.
5. The rest for an undecided link; doubling where the rule calls for it; then stage 4.
6. The summary, the finding, a ledger entry for the calls made running it, and the result to the maintainer
   with the collapsed-move plan's open calls.

If the move lands on bartcore before this runs, its builds are bartcore's tip and the tip flipped.

## Multiple comparisons

The maintainer placed the multiple-comparisons problem in a multiple-treatments scenario (2026-10-17: "It was
in a multiple treatments scenario."), though not the exact case, which their colleagues hold. Group M reads:
- A6, the maintainer's case: per unit, each null contrast's coverage of zero and the share of contrasts
  declared (M1); the share of replications in which any null contrast of any unit is declared, by the band
  over all units' contrasts and pointwise (M2); the band over one unit's six contrasts, its coverage reported;
  the share of declared contrasts with the wrong sign (M3). Its M1 and M2 alone are guards of their own.
- The null cells of A2 and A3: the share of units declared (M1) and the family-wise rate (M2).
- The sparse-heterogeneity cell: the same at its units with no effect, where false discovery sits beside real
  effects.
- Among declared effects in every non-null cell, the wrong-sign share (M3).

The colleagues' exact case, if reported, can be scored on the stored draws or run as a new cell after the
study.

## Amendments before results

Made by the orchestrator while the harness was written, 2026-10-17, before any pilot number was read:
- **R at the pilot.** The 26 B references a link, 8 chains of each build as above, run before the pilot
  reading. They are inputs, not scored fits. (R's margin, frozen from them in the first form of this
  amendment, went with dec-B456.)
- **M2's cells.** The 7 null cells and all four A6 cells, 316 pairs. The sparse cell's family-wise rate
  is reported beside M2, outside the rule.
- **T's cells.** 23 non-null cells in A2 to A4 (A2 3, A3 12, A4 8), the count the planned 692 pairs
  imply; "15" above is a miscount.
- **The large datasets.** The four of more than 5,000 rows (adult, bank, magic, mushroom), the k study's
  rule; "the large five" above is a miscount.
- **Order.** The A1 truths and the pilot's propensities run before the timing pass, so A1's timing cells
  draw their outcomes from the truths and the effect cells carry their real propensity.
- **The pilot** is replications 1 to 10 of every cell (B seeds 1 to 10, adult's 8), reused in the
  main run. Its standard error is the per-replication spread read from it, scaled to the planned pair
  counts, as the planned column was computed.
- **R is a diagnosis, outside the rule** (dec-B456): neither a gain group nor a guard, with no margin. Its
  numbers are reported under diagnosis, as whether a default-length chain reaches its own long run's
  posterior. The gain groups are F, T and M. The B references still run, before the pilot, as that
  diagnosis's input. The halfway reading waits on the A cells alone. Doubling is skipped when no gain
  measure's upper bound reaches its bar. The power statement drops R. The guards' RMSE and interval-score
  ratios pool the A cells. The equal-time choice pools the cells of F, T and M1.
- **Effect fits** call bart() with bartCause's default response model (the outcome on the covariates,
  the treatment and the estimated propensity, 10 chains, every unit predicted under both treatments),
  not bartc().

Made after the pilot was read, 2026-10-17, on the pilot's review (scratch/bis-pilot-rev-report.md). The
halfway look's F, A1's cells 1 to 8, is pilot data, so every change below was made with that F in view:
- **R on interior rows (fix).** The references and the stored draws were single-precision probabilities,
  which round a probability near 0 or 1 to exactly 0 or 1, and the mean over rows dropped every row where
  one column was undefined. R is now read on the rows or units whose reference interval (5th to 95th
  percentile) lies in [0.01, 0.99], an effect when both of its arms' probabilities do; the saturated rows'
  reference mass, and their share, are reported apart; each column is averaged over the rows where it is
  defined; R is reported with and without the flagged problems. References made from now on (stage 4) and
  the B fits' stored long draws keep the latent scale, every test column, and form the scored quantity in
  double precision when read, A6's as its six contrasts. The pilot's B references hold the scored quantity:
  there a probability is judged by its own interval, and an effect, whose arms it cannot see, is interior
  when its reference interval has positive width.
- **Near-deterministic cells apart (protocol, the maintainer rules).** Each A1 cell and link is classified
  from its truth alone, over every row of its base data: near-deterministic, the plan's own type, when the
  truth lies in [0.01, 0.99] on fewer than half the rows. Every reading reports these cells apart as a
  diagnosis, with the coverage, interior-row coverage, interval score and RMSE of both builds. Whether F
  and the A1 cells of the simultaneous guard leave them out is one switch in the rule's reading, off by
  default, as frozen; the pilot is read both ways, with the classification's sensitivity to the row cutoff
  (0.005, 0.01, 0.02) and the cell threshold (25, 50 and 75 percent interior). Nothing else in the rule
  changes.
- **The pilot's reading (fix).** The pilot's replications can be read as the halfway look would read them,
  for the record; the pilot spends no alpha.

## Calls made in planning

The orchestrator's, on the critique's open items, 2026-10-17:
- Equal time is measured from sampling time, `n.threads = 1` on every fit, the draws-or-burn-in choice made on
  the pilot and frozen, the simultaneous measures thinned to the move's draw count.
- The A1 BART truth is made by the move build at reference length; the additive truth is added.
- Twins and ACIC 2019's binary practice datasets are fetched; the high-dimensional, sparse, 1 percent,
  rare-under-weak-overlap and full-size adult cells are added; the multiple-treatments design A6 is added as a
  main problem in group M, the maintainer's case.
- Effect fits use 10 chains with one propensity per replication shared by both arms.
- The study runs under the bench-box grant with the stop raised to cover the doubled path.
- No rebase; stage 1 checks the study's paths on x86.
- Leave out only on harm shown; margins at 3.3 standard errors or more; gain group M added; a half still
  inconclusive after doubling stays out by default; references of 8 chains of each build, failing problems
  kept and flagged.

The maintainer's: the bar, dec-B455; the multiple-comparisons scenario, multiple treatments (2026-10-17).

Made here, beyond the critique:
- Replications cut to pay for the added cells within the stop: A1 from 20 to 16 a cell, A2 and A3 at 1,000
  units from 40 to 30, A3 at 4,802 units from 20 to 16, B1 from 20 seeds to 12 (adult 8) and B2 to 15, the cross-link reading from 5 to 3, diagnosis
  references from 3 replications a cell to 1 (2 in thin types), pooled by type. The counts above still
  resolve 0.02.
- References are not continued past their planned length: a failing problem is flagged instead.
- Doubling is skipped where it cannot change the outcome.

Declined from the critique:
- Rebasing first (its build point): the orchestrator ruled it out, and the drift probe found the shipped build
  equal to bartcore's tip.
- A check at 2 threads: every fit runs one thread, so threads are on no path the study takes.
- Using the generated true propensity: the estimated one is what a bartCause user's fit carries.
- Pooled signed measures in place of sized margins (offered as an alternative): the margins answer the power
  problem without changing what is measured.

## Probes run in planning

On the laptop under a load of 15 to 50, private library scratch/libs/bis (the branch tip, aciccomp2016 and
aciccomp2017), one job at a time:
- the build of 207f740e installs; the ACIC generators install offline from the local checkout;
- data: ACIC 2016 (4,802 by 58, 858 treated in its first setting's first replication, the generator 0.25 s);
  IHDP (747 after the deletion, 139 treated); lindner (996, 26 deaths); lalonde (614, 77 percent with
  earnings); ACTG 175 (2,139, 4 arms); the 16 UCI datasets from the cache; Twins (71,345 pairs, 54 columns;
  11,984 pairs under 2,000 g); ACIC 2019's practice sets (8 of 16 binary, sizes as above);
- one default fit under each link on IHDP with a binary outcome (1.7 s and 2.1 s);
- a 2,000-draw run's first 500 draws equal a 500-draw run's bit for bit under both links, and threads 1 and 2
  draw alike;
- a sampler from bart() with `samplerOnly = TRUE` does not draw bart()'s chain at the same seed;
- duplicate covariate patterns in the UCI datasets.

The critique's probe (scratch/bisc-crit.md): the shipped build against bartcore 92ae2714, bit for bit on
arm64.

## Not verified

- The shipped build against bartcore 92ae2714 on x86 (stage 1).
- The x86 time of the effect designs, the 10-chain fits, the high-dimensional cells, Twins and adult at full
  size; the cost table scales the k study's and csmx's timings.
- The per-replication spreads the margins are sized on (the pilot rereads them).
- That references reach acceptance, probit near separation least of all, and how many problems are flagged.
- That the additive truth fits every dataset (tictactoe and mushroom are all factors and nearly separable).
- The memory peak of the largest jobs on the box.
