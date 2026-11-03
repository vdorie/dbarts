# factor-split-probability: should a factor's prior split probability grow with its number of levels?

Status: CLOSED 2026-11-03 without a run (dec-B576): a factor keeps one column's share of the split prior,
since a prior should treat variables alike and no source argues a factor deserves more. Planned and amended
twice the same day; kept for its survey of what the literature and other packages do and its probes. Nothing
ran beyond the sizing probes below.

agent: the blind critique (opus) ran; this revision answers it. One opus implementer writes the harness,
runs the smoke on the laptop and the main run on the bench box, and drafts the report; one opus review of
the harness and smoke before the main run, since the harness carries the statistics; two independent
referees (opus) on the report before it reaches the maintainer.
rng: NEUTRAL. No package code changes: every candidate rule is a `split.probs` vector the shipped package
already takes, and the BayesTree reference arm is `factors = "indicators"`. Only the conditional DART stage
needs a scratch build (about 30 lines, never landed). A change of default the study supports is a separate,
draw-changing slice with its own plan.
window: after the merge to main, before 1.0-0 (dec-B78, VD 2026-09-28).
budget: harness about 800 lines of R (benchmarks/R/factor-split-probability.R) reusing the house's binary
dataset loaders. Compute about 9 box core-hours if no candidate passes the real-data reading, 20 to 46 if
one does (the treatment-effect guard is sized to its margins), stop at 65 ([Sizing](#sizing)); report about
5 pages of text plus figures.

## Summary for the maintainer

**The question.** When a tree picks a variable to split on, should a factor with many levels be picked
more often than a numeric column, by some function of its number of levels?

**Why it is open.** The BayesTree-style coding turns a factor with K levels into K indicator columns, each
picked as often as a numeric column. Since 1.0-0 a factor is one column split on any subset of its levels
(dec-B78), picked as often as any other column. You asked whether the weighting the old coding gave should
come back on the new one. "The weighting it gave" has two readings, and they differ a lot: the old coding
gave a factor K times a numeric column's share of splits, but each of those splits only separates one
level from the rest. Counted by how often a split separates two rows at random, the old coding gave a
factor 4 to 6 times a numeric column's power, the new coding gives 1.3 to 1.5 times, and putting K's
weight on the new coding would give 4 times at 3 levels and 74 times at 50.

**What others do.** Every package that sets a factor's weight on purpose gives it one column's weight: the
original Bayesian CART prior (pick a variable uniformly), bartMachine (each indicator weighted 1/K, since
its version 1.0.2), stochtree, SoftBart's formula interface, Linero's grouping prior, flexBART and
pymc-bart. Only BayesTree and the BART package weight by level count, both only because they code factors
as indicators; the BART package computes a grouping of a factor's indicators that its sampler never reads.

**What the probes found.** Weighting by level count cut the error of the fitted function by 0 to 8 percent
where one 10-level factor drives the function (the four published flexBART designs), raised it by about 50
percent where factors carry no signal (nearly a third of the splits went to a useless 50-level factor), and
made no clear difference on three small real datasets. With no signal, each predictor's share of splits
after fitting equals its prior share, so any weighting biases split-count importance toward many-level
factors by exactly its weight. Which rule is better therefore depends on the data users meet, and the
study reads that on real data.

**Candidates.** Today's rule (one column); the level count (K, 1 for a binary factor); a weight matched to
the old coding's separating power (3 at 3 levels, rising to 4); and the square root of K - 1, the middle
of the range between them. The old indicator coding runs beside them as a reference only.

**Data.** Fourteen of the sixteen real regression datasets the flexBART paper assembled for factors in
BART (2 to 84 levels; two left out because their response is almost exactly a function of the predictors,
named below); seven real binary datasets with factors from the house's set; then, only if a candidate
passes on real data, the 2016 causal inference competition's data, real covariates with known treatment
effects. Small published simulations explain the mechanism and do not decide.

**What decides, fixed before any fit.** Today's rule stays unless a candidate is practically better across
the real datasets: a gain of more than 1 percent in held-out error, the margin and the test (a Bayesian
signed-rank test across datasets) both taken from a standard tutorial on comparing methods over many
datasets (Benavoli, Corani, Demsar and Zaffalon 2017), with no single dataset made practically worse by
the same tutorial's test within a dataset; a loss on any one dataset makes the result a trade-off, and
today's rule stays. The margin is on the error of predicted outcomes, not of the fitted function. A
candidate that passes must then lose no more than 1 percent in error, and no more than 2 points of
interval calibration (your bar from the binary interval study, over- and under-coverage both counting),
on the competition's treatment effects; their bias is reported. A result of no difference is reported as
such, not only as a failure to find one.

**Your call** (below): whether the bias in split-count importance that a weighted rule brings by
construction is a reason on its own to keep today's rule. Recommended: no; report it, and say it in the
manual if a weighted rule wins.

**Cost.** About 9 core-hours on the box if no candidate passes on real data; 20 to 46 if one does, most of
it the treatment-effect check, sized so that a candidate as good as today's rule can pass it; stop at 65.
The result comes to you as a short PDF report after two independent referees.

## Goal

A research report (a PDF of at most five pages of text plus figures, self contained) that answers whether
the default `split.probs` for a single-column factor should depend on its level count, argued as the BART
and method-comparison literature argues, refereed twice, with a recommendation the maintainer can rule on.
Nothing in the package changes under this plan.

## Context

- dec-B78: "Weighting a factor column's split mass by its number of levels is a research item, moved on
  2026-09-28 from after the release to between the merge and the submission ("I would actually move up
  the research on split probability to post-merge, pre-release"). The maintainer chose "Subset splits as
  the modern default."" The encoding is ruled; only the weighting is studied.
- dec-B559: the study "hold[s] up to academic scrutiny", comes as "small research reports (pdf, at most 5
  pages of text with as many graphs as needed, self contained and easily read)", and imitates "how other
  people argue for default priors, not inventing things ourselves." docs/plans/README.md, Studies, adds the
  dec-B558 questions: datasets people meet, natural estimands, RMSE, bias and heterogeneous effects.
- dec-B512 (revising dec-B116, dec-B246, dec-B252): "split.probs takes a vector of split probabilities or
  dart(), which is DART, a prior on which variable a split uses". A changed default is the value
  `split.probs` takes when unstated, at each door that takes it.
- dec-B338: "split.probs weights how often each predictor is picked for a split"; negative, non-finite and
  all-zero vectors are refused. dec-B547: a named vector reads by column or term name with an optional
  `.default`.
- dec-A09, dec-A115: an ordered factor splits between its level codes, as a numeric column does, from its
  declared levels. Out of scope: it is weighted as a numeric column is.
- dec-B370, dec-A26: wide factors are stored sparse and split by subset like any factor; the rule applies
  to them unchanged.
- The engine: `CGMTreePrior::drawSplitVariable` picks a column among those available at a node with
  probability proportional to its `split.probs` entry (uniform when none), and
  `CGMTreePrior::drawRuleForVariable` then draws a factor's rule uniformly over the left-right patterns of
  its reachable levels that leave neither side empty. `DartPrior` draws the probabilities from a symmetric
  Dirichlet, alpha / p per column plus the split counts.
- The help (man/bart.Rd, `factors`) already says subset splits "give a factor the prior split probability
  of one predictor rather than one per level".

### What the literature and other packages do (read from their sources)

| source | encoding | a K-level factor's share of the split prior |
|---|---|---|
| Chipman, George, McCulloch (1998, 2010) | as given | "the uniform prior on available variables" (CGM 2010, sec 2.2, following CGM 1998); the 2010 bake-off converted categoricals "into 0/1 indicator variables corresponding to each level" (sec 5.1) |
| BayesTree 0.3-1.5 (`makeind`) | K indicators, 1 for K = 2 | K, 1 for a binary factor |
| BART 2.9.10 (`bartModelMatrix`, `class.ind`) | K indicators, also for K = 2 | K; `grp` is computed and passed to C++, which never reads it |
| dbarts 0.9-34; in 1.0-0 `bartBT` and `factors = "indicators"`, which keep its expansion (dec-A09) | K indicators, 1 for K = 2 (probed on 1.0-0) | K, 1 for a binary factor |
| bartMachine 1.4.2 (`cov_prior_vec`, since 1.0.2) | K indicators | 1: each indicator weighted 1/K ("downweighting factors by their number of levels") |
| stochtree 0.4.5 (`variable_weights`) | one-hot | 1: the original variable's weight divided among its indicators |
| SoftBart 1.0.3 (formula path, `group`) | K indicators | 1: one Dirichlet share per factor, uniform among its indicators |
| Linero and Yang (2018), grouping prior | indicators | 1: u over groups, v within; "an extension of the approach used by the bartMachine package" |
| flexBART 2.0.6 (Deshpande 2024) | one column, each level to a side with probability 1/2 | 1, under a sparse Dirichlet by default |
| pymc-bart (`split_prior`, `SubsetSplitRule`) | one column | 1 |
| dbarts 1.0-0 | one column, subset splits | 1 |

Variable selection bias in trees. CART's exhaustive search favors predictors with more candidate splits
(Breiman, Friedman, Olshen and Stone 1984); unbiased selection was the aim of Loh and Shih (1997), Kim and
Loh (2001) and Hothorn, Hornik and Zeileis (2006); Strobl, Boulesteix, Zeileis and Hothorn (2007) show the
same bias in random-forest importance, holding that "in the null case ... a sensible variable importance
measure should not prefer any one predictor variable over any other", with X1 ~ N(0, 1) and factors of 2,
4, 10 and 20 equally likely categories at n = 120, a null case and a power case in which a binary response
depends on X2 with relevance 0.05, 0.1, 0.15 or 0.2. Their mechanism is search; here the bias would come
from the prior weight itself, a different mechanism with the same symptom (probe p5). In BART, Bleich,
Kapelner, George and Jensen (2014) read variable inclusion proportions against a permutation null.
Deshpande (2024) shows one-hot coding cannot pool levels and that uniform subset splits are "biased
towards balanced partitions and ... generally unable to isolate a single outlying level" (his DGP2); his
benchmark compares flexBART's coding and weight with BART's coding and weight together, so it cannot say
which of the two did the work. This study separates them.

Comparing methods over many datasets. CGM 2010 (sec 5.1) argued its defaults by out-of-sample RMSE over a
benchmark of real datasets, 20 random splits each, summarized as relative RMSE. Deshpande (2024, sec 5.1)
used one-sided paired t-tests per dataset with a Bonferroni correction over 32 comparisons. Benavoli,
Corani, Demsar and Zaffalon (2017, JMLR) replace such tests with Bayesian ones over datasets with a region
of practical equivalence: "it is sensible to define that two classifiers whose mean difference of
accuracies is less that 1% are practically equivalent", a Bayesian signed-rank test (Dirichlet process,
prior strength s = 0.5, pseudo-observation at 0), and a decision when the posterior probability of a
region exceeds 0.95; their correlated t-test is Nadeau and Bengio's correction for overlapping training
sets. See also Benavoli, Corani and Mangili (2016) on post-hoc tests over mean ranks.

## Probes (laptop, shipped build, 2 cores; scratch/fsp-plan/probes, scratch/fsp-crit/probes)

| probe | design | finding |
|---|---|---|
| p1 | flexBART DGP 1-4, n = 1000, 10 reps, 2 chains | log RMSE ratio, level count vs uniform: -0.07, -0.04, -0.08, 0.00 (sd per rep 0.06 to 0.07); indicators vs uniform +0.11, -0.09, +0.11, +0.06 (Deshpande's pattern); coverage of f, paired difference -0.01 (sd 0.02); posterior split share of the factor 0.24 uniform, 0.51 weighted, against prior 0.09 and 0.50 |
| p2 | CGM's Friedman function on 5 of 10 numerics plus null factors of 2, 4, 10, 20, 50 levels, n = 500, sigma 1, 10 reps | log RMSE ratio +0.41 (sd 0.13) for level count, +0.33 indicators; coverage 0.934 uniform, 0.898 level count; the null 50-level factor took 0.30 of the splits under the level count, 0.026 under uniform |
| p3 | real: ggplot2 mpg, carData Salaries, MASS Cars93; 20 splits of 75/25 | log RMSE ratio level count vs uniform +0.04, -0.003, +0.05, per-split sd 0.20, 0.007, 0.16 |
| p4 | ACIC 2016 covariates, a toy surface, default `bart` | 5.6 s per fit on 2 threads (about 11 core-seconds); level count and uniform within 5 percent in time |
| p5 | Strobl null, gaussian response, n = 120, 20 reps | posterior split share equals the prior share: uniform 0.20 each; level count 0.03, 0.03, 0.11, 0.28, 0.55 against prior 0.03, 0.03, 0.11, 0.28, 0.56; sd per rep under 0.008 |
| separation (critic) | prior probability that one split separates two random rows, equal level frequencies, relative to a numeric column | subset split 1.50, 1.33, 1.29, 1.28, 1.35, 1.43, 1.47 at K = 2, 3, 4, 5, 10, 20, 50; indicators 1.50, 4.00, 4.50, 4.80, 5.40, 5.70, 5.88; level count on subsets 1.50, 4.0, 5.1, 6.4, 13.5, 28.5, 73.5 |
| power (critic) | the first plan's R1 (Wilcoxon over 16 datasets, Holm over 3) on simulated per-dataset means | a uniform 0.5 percent gain passes 15 percent of the time; 5 percent on 6 datasets, 23 to 26 percent |
| bsr | the amended R1 (Bayesian signed-rank, rope 0.01) over 14 datasets at the critic's per-dataset SEs, 400 simulations | declares a gain: 0.86 at a uniform 2 percent gain, 0 at 0.5 percent, 0.38 at 5 percent on 6 of 14, 0.87 at 5 percent on 8 of 14, 0 under no difference; declares equivalence 0.57 of the time under no difference |
| timing (critic) | n = 18,000, factors of 84, 33, 18, 6 and 5 levels | 15 to 17 core-seconds a fit, A and B alike |

## Candidates

Weights apply to an unordered factor's column when `split.probs` is unstated; numeric, logical and
ordered-factor columns keep weight 1; K is the number of levels with training rows; every rule gives a
binary factor weight 1.

- A, uniform (today): w = 1. Split counts compare fairly across columns under no signal (p5). A
  many-level factor that matters is proposed 1/p of the time.
- B, level count: w = K for K >= 3. The old coding's share of split proposals moved onto subset splits. It
  is not the old coding's effect: each subset split separates more than an indicator split does, so B
  gives a factor 4 times a numeric column's separating power at K = 3 and 74 times at K = 50, where the
  old coding gave 4 to 6 times (probe separation). A 50-level factor beside 10 numeric columns takes 83
  percent of the prior splits.
- M, matched to the old coding's separating power: w = 4 (1 - 2^(1 - K)) for K >= 3, that is 3, 3.5,
  3.75 at K = 3, 4, 5 and about 4 beyond. It is the ratio of the indicator coding's pair-separating power,
  2 (1 - sum p_k^2), to the subset split's, (1 - sum p_k^2) 2^(K-2) / (2^(K-1) - 1), so it holds for any
  level frequencies p_k: the frequency term cancels. It is the literal reading of what the old coding
  "gave implicitly" on the new one, and the only rule here that does not grow without bound in K. It was
  derived in this study, not taken from a source, and the report says so.
- S, sqrt(K - 1). No source argues it. It lies between A and B for every K and crosses M near K = 17, so
  a gain at B and a loss at A cannot both be missed through an unsearched middle of the family.
- R, reference only: `factors = "indicators"`, the old coding and weight together, on the real regression
  block and the competition data's 10-replication pilot. Not a candidate (dec-B78 rules the coding); it says
  what 0.9-34 users had.
- Dropped after the critique: C, w = K - 1, which differs from B only at K = 3 and 4.
- Considered and not fitted: the log of the number of distinct subset splits, log2(2^(K-1) - 1), within
  0.2 of K - 1 for K >= 4; a node-local weight from the levels reaching the node, an engine change reached
  only if B wins; and BART's K = 2 weight of 2, an artifact of its coding that its documentation
  contradicts ("one dummy is created if q=2").

**With a stated `split.probs`.** A stated vector is literal (dec-B338) and no rule applies. How a named
vector's `.default` reads under a weighted default is an implementation-time call.

**With DART.** `dart()` takes no vector; its Dirichlet is symmetric, alpha / p per column. A weighted
default carries to DART as a non-uniform base measure, alpha * w_j / sum(w) per column. Since the split
counts dominate the Dirichlet once trees use a column, the weight is expected to matter less there. Read
only if a candidate passes R1 and R2 (step 7), on a scratch build of that base measure.

**With forest column restrictions and interaction constraints.** The weights renormalize over the columns
a forest or node may use, as today.

**Mixing against prior.** A weight changes the prior and how often the sampler proposes the factor at
once; at equal running time the study measures both, as users get them. For a candidate that passes R1,
long chains say which a gain comes from; a gain from mixing alone would point to a proposal change that
keeps the uniform prior, a separate design, not to a new default.

## Estimands and scores

Every summary is computed per draw and averaged last. Intervals are central 95 percent posterior
intervals for f and effects and 90 percent predictive intervals for new outcomes, scored as distance from
nominal, over- and under-coverage both counting; width is reported beside coverage.

- f(x) where the truth is known: RMSE; bias as defined for simulation studies (Morris, White and Crowther
  2019): at each fixed test point, the mean over replications of the posterior mean minus the truth,
  reported as its root mean square over test points, so that errors of opposite sign at different points
  do not cancel; coverage and width.
- Level means, the average of f over the rows at each level, the quantity an analyst reports for a
  factor: RMSE, root-mean-square bias over levels and coverage.
- Treatment effects (ACIC 2016), as Dorie, Hill, Shalit, Scott and Cervone (2019) scored the competition:
  the average effect on the treated (bias, reported; RMSE, coverage and width), individual effects (PEHE and
  the coverage of their intervals), and the effect on the treated within each level of the 16-level factor
  x_21, each level weighted by its treated rows. x_21's levels hold 8 to 2,561 rows, so unweighted levels
  would let a level of 8 rows count as much as one of 2,561; weighting by rows scores the effect per treated
  individual, the competition's own unit.
- Real outcomes: held-out RMSE (CGM 2010; Deshpande's SMSE ratio is its square), summarized also as CGM's
  relative RMSE; coverage of predictive intervals; for binary outcomes, held-out log loss and Brier score.
- Variable importance: each predictor's share of splits (CGM 2010 sec 5.2 reads usage on the Friedman
  data) and the share of draws using it; Strobl's power case, the informative predictor's share and rank;
  on ACIC, the treatment column's share of splits (its prior share falls from 1/59 to 1/83 under B).
- Time per fit.

**Equal time.** Fits run at the shipped defaults (4 chains, 500 burn-in, 500 draws, 75 trees). A candidate
that passes R1 and is slower than A at the median time per fit has its comparison re-read with A rerun at
burn-in and draws scaled by the time ratio.

## Data

1. **Real regression datasets with factors**: the sixteen of Deshpande (2024), Table 1, with his sources
   (UCI, the Journal of Applied Econometrics data archive, R packages, fueleconomy.gov, StatLib, OzDASL).
   Disclosed: they were assembled by a method developer to show his method's handling of factors, so they
   lean toward data where factors matter, the population the question is about but not every user's data.
   Two are left out before any fit, by named reason, not by a number: Caschool, whose response (testscr)
   is the mean of two of its predictors (readscr, mathscr), and fuelEcon, on which all three of
   Deshpande's methods reach an SMSE of 0.000; on a response that is a near-exact function of the
   predictors, held-out error measures how finely trees approximate a deterministic surface, not how a
   prior trades signal against noise. That leaves fourteen (abalone, ais, Alcohol, amenity, attend, cane,
   cpu, Engel, Insur, Medicare, mpg, servo, spouse, strike). Step 1 fetches each, caches it with a recorded
   sha256 as benchmarks/R/uci-binary.R does, and checks dimensions and level counts against Table 1; one
   that cannot be fetched or matched is recorded and left out, with no substitute. The block decides only
   if what remains still holds both a dataset with a factor of 20 or more levels and one whose factors all
   have fewer than 6, since the question is about level counts and a block without both ends cannot speak
   to it; otherwise it is reported and does not decide. 20 random 75/25 splits each (CGM 2010's count;
   under the overlap correction more splits barely shrink a dataset's error), test rows at a level absent
   from training dropped and counted (dec-A85), the same rows for every arm.
2. **Real binary datasets with factors**: the house's binary datasets (benchmarks/R/binary-datasets.R,
   uci-binary.R; assembled for the binary k studies before this question) that have an unordered factor
   of 3 or more levels, less two left out by named reason: tictactoe, whose nine predictors are all
   3-level factors, so every rule gives every column the same weight and the arms are the same model;
   and mushroom, whose classes are separated almost exactly by single predictors, so held-out log loss is
   near zero for every arm. That leaves adult, bank, german, creditapproval, cleveland, infert and birthwt
   (step 1 confirms each has such a factor). The span condition of Data 1 applies here too, since this
   block can pass R1 on its own. The house's split and cap: 80/20, at most 4,000 training rows past 5,000
   rows. 20 splits each.
3. **ACIC 2016**, run only for a candidate that passes R1 (Dorie et al. 2019; aciccomp2016): 4,802 rows of
   real covariates, 58 columns with factors of 6, 16 and 5 levels (x_2 and x_24 each have one level holding
   98 and 72 percent of the rows); all 77 published settings, taken as fixed: they are the competition's
   whole design, a grid of the organizers' settings (Dorie et al. 2019) rather than a sample from a
   population of settings, and the competition reported its scores over them. So the estimand is the average
   over these 77, and its error comes from replication within settings alone; the competition's potential
   outcomes are the truth. A fit puts the treatment as a predictor and reads effects from the two
   counterfactual predictions (Hill 2011). One covariate set: a guard against harm to effects, not evidence
   that a gain generalizes.
4. **Designed, to explain the mechanism (they do not decide)**: flexBART DGP 1-4 (Deshpande 2024, sec 4.1;
   ten uniform numerics and one 10-level factor driving the function), n = 1000, 25 replications, 500 test
   points fixed across replications; Strobl et al. (2007), the power case at its four relevance values and
   a 20-replication null (p5 already shows the null's answer), binary response, 50 replications; CGM
   2010's Friedman function on 5 of 10 numerics plus null factors of 2, 4, 10, 20 and 50 levels, n = 500,
   25 replications, test points fixed; for a candidate that passes R1, long chains on the flexBART DGPs,
   A and that candidate, 10 replications at four times the default length.

## The decision rule (fixed before any fit)

For candidate c (B, M or S) and real dataset j, L_cj is the mean over splits of log(RMSE_c / RMSE_A); for
binary data, of the ratio of held-out log loss. The region of practical equivalence is [-0.01, 0.01], a
1 percent relative difference in error, Benavoli et al.'s 1 percent carried from accuracy to relative
error. It is a margin on the error of predicted outcomes, which on noise-dominated data (Alcohol's SMSE
is 0.96) is a far larger relative change in the error of the fitted function; the report words
"equivalent" as equivalent in predicted outcomes, not in f.

- **R1, real outcomes**, by Benavoli et al.'s Bayesian signed-rank test (s = 0.5, pseudo-observation at 0)
  on the L_cj of each block: c passes if in at least one block the posterior probability that c is
  practically better (the region where that probability is the largest of the three) exceeds 0.95, and in
  neither block does the probability that c is practically worse exceed 0.95; and no single dataset is
  practically worse: by the same tutorial's correlated Bayesian t-test within each dataset (posterior
  Student t on the split-level log ratios, n - 1 degrees of freedom, scale sqrt(1/n + rho / (1 - rho))
  times their sd, rho the test share of the rows), any dataset whose probability of L_cj above 0.01
  exceeds 0.95 makes the result a reported trade-off, not a pass. Among several that pass, the one with
  the larger probability in the regression block. No multiplicity correction: per candidate, a false pass
  under no difference ran 0 at precise per-dataset errors and 1.5 percent (regression) and 0.7 percent
  (binary) at noisy ones (probe bsr2), so the three candidates' family-wise rate is at most about 7
  percent at the noisy errors by the union bound, near the conventional 5. CGM's relative-RMSE boxplots
  are reported, not decisive. A candidate whose probability of practical equivalence exceeds 0.95
  is reported as equivalent to today's rule; otherwise as inconclusive.
- **R2, treatment effects (ACIC 2016)**, for a candidate that passed R1: non-inferiority, measure by
  measure, over the 77 fixed settings. For each measure the difference c minus A is averaged within each
  setting, then over settings with equal weight; its standard error is from the within-setting variances
  alone, sqrt(sum_s sigma_s^2 / R) / 77. c passes a measure if the one-sided 97.5 percent upper bound
  (the upper end of a two-sided 95 percent interval, the usual non-inferiority convention, ICH E9) lies
  below the measure's margin:
  - error measures, the log ratio of SATT RMSE, of PEHE and of the x_21 per-level RMSE: 0.01, the margin
    of R1 (Benavoli et al.). SATT's bias is inside its RMSE and is reported, not guarded separately;
  - calibration measures, the coverage of SATT, of individual effects and of f(x), each scored as a
    setting's distance from 0.95 corrected for the noise in its coverage as the binary interval study
    corrects it (dec-B455): 0.02, the maintainer's bar for a coverage gap in that study (dec-B455). A
    1 percent margin carried from error to coverage points had no argument, and on SATT's one yes-or-no
    outcome per fit it would need about 1,000 d replications per setting, d the rate at which the arms
    disagree, which the budget cannot reach.

  c passes R2 only if it passes every measure, an intersection-union test, which needs no multiplicity
  correction (Berger 1982). Replications: the first 10 per setting are an internal pilot that re-estimates
  only the within-setting variances (Wittes and Brittain 1990); the count R for A and c is then the
  smallest giving 80 percent power at a true difference of 0 on every measure, R = 7.85 sum_s sigma_s^2 /
  (77 m)^2 for margin m, capped at 40 by budget; a measure short of its power at the cap leaves R2
  inconclusive.
- **Verdict.** A candidate passing R1 and R2 is recommended. Otherwise today's rule stays: if no candidate
  passes R1, reported with each candidate's equivalence reading; if one passes R1 and fails or cannot
  settle R2, reported as a trade-off.
- The designed cells, variable importance and the long chains are reported beside the verdict to explain
  it. Whether the importance bias of a weighted rule is decisive on its own is the open call.

## Sizing

What R1 detects, from probe bsr2 (the signed-rank test with the per-dataset harm guard, 400
simulations) at two per-dataset precisions: precise, SE 0.002 to 0.05, and noisy, SE 0.004 to 0.12 (the
critic's re-check; p3's mpg-like data reach 0.12). A dataset's SE under the overlap correction is 0.62
times its per-split sd at 20 splits and 0.59 at 50.

| component | question | size | box core-h | what it resolves | sequential stop |
|---|---|---|---|---|---|
| real regression | R1 | 14 x 20 splits x 5 arms (A, B, M, S, R) = 1,400 fits | 6.6 | passes, precise / noisy: uniform 2 percent gain 0.87 / 0.39; uniform 5 percent 0.99 / 0.92; 5 percent on 8 of 14, 0.82 / 0.55; on 6 of 14, 0.37 / 0.38; no difference 0 / 0.015; one dataset losing 49 percent beside 13 gaining 2, 0 / 0 (the guard fires every time, where the signed-rank test alone declared a gain 0.62) | none: more splits do not move a dataset's corrected SE |
| real binary | R1 | 7 x 20 x 4 arms = 560 | 1.1 | passes, precise / noisy: uniform 2 percent 0.42 / 0.24; 5 percent 0.95 / 0.69; no difference 0 / 0.007 | none |
| flexBART DGPs | explain | 4 x 25 x 4 arms = 400 | 0.4 | log RMSE ratio SE 0.014 (p1 sd 0.07): 3 percent; paired coverage SE 0.004 | stop at 25 |
| Strobl | importance | 4 relevances x 50 + 20 null, x 4 arms = 880 | 0.5 | split share SE about 0.001 (p5 sd 0.006); the informative predictor's rank | stop at 50 |
| Friedman + null factors | explain | 25 x 4 arms = 100 | 0.1 | log RMSE ratio SE 0.03 on an effect of 0.4 | stop at 25 |
| ACIC 2016 (on an R1 pass) | R2 | pilot 77 x 10 x 3 arms (A, c, R) = 2,310; then A and c to R reps, at most 77 x 30 x 2 = 4,620 more | 10.6 to 31.8 | R = 1,019 sigma^2 for an error measure (per-fit sd sigma of the log ratio: 10 reps at 0.1, 41 at 0.2) and 255 d for a coverage measure (13 reps at d = 0.05, 38 at 0.15) | the count set once from the pilot |
| long chains (on an R1 pass) | mixing or prior | 4 x 10 x 2 arms, 4x length | 0.4 | the p1 gap of 0.07 at SE 0.02 | none |
| DART (on an R1 and R2 pass) | does the weighted base measure hold | A and c under DART, real blocks and 10 ACIC reps | about 5 | as above | none |

Box core-hours are laptop core-seconds times 1.5, a ratio not yet measured that the smoke re-prices; fits
at the default 4 chains. About 9 if no candidate passes R1 (8.7, run first); 20 to 41 if one passes, 25
to 46 with the DART stage, near the top if PEHE's per-fit spread is near 0.2; stop at 65.

## Steps

1. Harness benchmarks/R/factor-split-probability.R: loaders for the fourteen regression datasets (fetch,
   sha256, dimension check against Deshpande's Table 1), the seven binary datasets from the house's
   loaders, ACIC 2016 through aciccomp2016, the designed generators with test points fixed across
   replications, the arms as `split.probs` vectors from each training set's level counts, one block per
   invocation writing one rds, and a summarize mode computing R1, R2 and the reported tables. Verify: every
   loader's dimensions printed against the table; two prior figures with no fits, every factor's share of
   the split prior and its pair-separating power relative to a numeric column under every arm for every
   dataset; the block-span condition of Data 1 checked.
2. Smoke on the laptop at 2 cores: one split or replication per block per arm; time per fit, re-pricing
   the table above.
3. Opus review of the harness and smoke (pairing, per-draw summaries, fixed test points, the Bayesian
   signed-rank and correlated t-tests against Benavoli et al.'s worked examples, the R2 bounds and pilot,
   the truth from ACIC). Fix round.
4. Pre-registration: commit the plan, any amendment from step 3, and the harness at one hash before the
   main run; later amendments are labeled as such.
5. Main run on the box, real blocks first, then the designed cells. Read R1.
6. For a candidate passing R1: ACIC replications 1 to 10 for A, c and R, the replication count from the
   pilot, the rest for A and c only, then R2; the long chains; the equal-time re-read if c is slower.
7. For a candidate passing R1 and R2: the DART stage on a scratch build whose base measure is
   alpha * w_j / sum(w); the node-local variant is named as follow-up work.
8. Report: R Markdown source docs/design/factor-split-probability-report.Rmd rendering every figure from the
   stored outputs; at most five pages of text. Outline: (1) the question, and the two readings of what the
   old coding gave; (2) what the literature and packages do; (3) the candidates and what each means for a
   user, with the two prior figures; (4) data and estimands, how each was chosen and what was left out and
   why; (5) results: the posterior simplex of the signed-rank test per candidate and block (Benavoli et
   al.'s figure), CGM-style relative RMSE boxplots by dataset, the binary block, the designed cells (error,
   pointwise bias and coverage by arm), the Strobl power and null shares, and where run, ACIC's effect
   scores and the long chains; (6) recommendation and limits. A finding note
   docs/design/factor-split-probability.md records the verdict with a link to the report.
9. Referees: two independent opus referees, each on a neutral brief (scratch/fsp-peer-common.md, modeled
   on scratch/nbs-peer-common.md: referee it as a methods-journal submission, owe it no deference, verify
   before asserting, read-only access to the stored outputs, small simulations of their own allowed,
   verdict first, ranked findings, the smallest addition that would change their verdict priced in
   core-hours). The brief names no preferred outcome. The report is revised to answer them; the maintainer
   receives the report, both referee reports and the response.
10. Records: plan Status and a Landing note, the finding, TODO "factor encoding", the INDEX rows; the
    ruling goes to the ledger when the maintainer makes it.

## Constraints

- No package code changes; NEUTRAL. A recommended change is its own plan and slice (draws change).
- At most 2 cores on the laptop; the box per its own limits; box folders deleted at landing, study data
  archived to the laptop first.
- No reference chains as truth; real outcomes scored on held-out rows only.
- The decision rule and the exclusions above are frozen at step 4; any amendment after a fit is scored is
  labeled as such.
- Out of scope: the coding (dec-B78), ordered factors (dec-A09), `bartBT`, numeric columns' weights, the
  subset-split rule's own distribution over partitions (Deshpande's DGP2 point), a proposal change that
  keeps the prior.

## Verification

- Step 1: `Rscript benchmarks/R/factor-split-probability.R datasets` prints each dataset's rows and factor
  level counts beside Deshpande's Table 1, mismatches marked, and the block-span check; `prior` prints both
  prior tables.
- Step 2: `Rscript benchmarks/R/factor-split-probability.R smoke DIR` completes every block, and
  `summarize DIR` refuses to mix smoke and main outputs.
- Steps 5 and 6: `summarize DIR` prints R1 (each candidate's three posterior probabilities per block) and,
  where run, R2 (each measure's mean difference and upper bound), with the counts of fits per cell, no cell
  short.
- Lint and air on the harness; `tools/check-doc-freshness.R`.

## Open calls for the maintainer

1. **Does split-count importance decide on its own?** Background: with no signal, the share of splits each
   predictor takes after fitting equals the prior share (probe p5), so any weighted rule makes split counts
   favor many-level factors by exactly its weight. This is the prior's own stated belief that such a factor
   is more likely to be split on, not the flaw of greedy search the tree literature names (Strobl et al.
   2007), though the symptom is the same; and a split count can be read against a permutation null, as
   Bleich et al. (2014) read BART's inclusion proportions, which removes the prior's share from the reading.
   dbarts reports split counts, and users read them as importance. Options: (a) importance decides alongside
   the estimates: every weighted rule fails it by construction, so today's rule stands without the study,
   saving about 9 core-hours (the cost unless a candidate passes on real data); (b) the estimates decide (R1
   and R2), importance is reported, and a winning weighted rule ships with the manual saying split counts
   then favor many-level factors in proportion to their weights; no extra cost; (c) the estimates decide,
   and a winning rule must also not lose the informative predictor in the Strobl power case; no extra cost,
   a rule that may fail on a 120-row design. Recommended (b): estimates and intervals are what users act on
   (your stated priority, 2026-10-17, that f and treatment effects come first), and a user who wants even
   counts can state `split.probs`.

## Amendments after the blind critique

The critique (verdict: run after amendments) made ten amendments; each was checked before it was taken.

1. B is not 0.9-34's prior beyond K = 3. Checked by rerunning the critic's separation probe: the subset
   split separates two random rows 1.3 to 1.5 times as often as a numeric column, the indicator coding 4
   to 6 times, B on subsets 4 to 74 times. B is redescribed, M added as the weight that matches the old
   coding's effect, and a separating-power figure added to step 1 and the report.
2. Drop C, keep S, add M. Taken: C differs from B only at K = 3 and 4. S is kept as the middle of the
   family, not by CGM's grid, which was per-dataset tuning; M is kept because without it the study cannot
   answer the TODO's literal question, the weighting the old coding gave, which B overshoots twelve times at
   50 levels.
3. R1 could move the default on a negligible gain and missed concentrated ones. Checked by rerunning the
   critic's power probe (a uniform 0.5 percent gain passed 15 percent of the time). R1 is now Benavoli et
   al.'s Bayesian signed-rank test with their 1 percent region of practical equivalence, which never
   declared a 0.5 percent gain in 400 simulations and declares a uniform 2 percent gain 0.86 of the time.
4. Sizing as what the rule detects, keep as an equivalence bound. Taken: the sizing table now states what
   R1 resolves, and R1 reports equivalence where its posterior supports it.
5. R2 margins per measure, per setting over the 77 settings, no multiplicity relief. Taken: a
   non-inferiority bound at the same 1 percent margin per measure, over setting means, every measure
   required (intersection-union, no correction); Bonferroni on a guard ran the wrong way. Where I differ
   from the critic's simpler study: its 10 replications per setting cannot certify non-inferiority at that
   margin unless the differences are small (at 5 replications the critic measured detectable losses of 2.6
   percent in PEHE and 1.6 points of coverage), so the count is set by an internal pilot and capped at
   40, which is why the pass path costs up to 42 core-hours on ACIC, not about 7.
6. Blind exclusions and a floor. Taken with named reasons, no numeric cut: tictactoe (arms identical),
   mushroom (near-deterministic), Caschool (response defined from two predictors), fuelEcon (every
   published method at SMSE 0.000). The floor is a span condition argued from the question, not a count,
   since the Bayesian signed-rank test can decide on as few as two datasets (probe bsr) and so supplies
   none. The method-developer origin of the sixteen is disclosed.
7. x_21 weighting and the treatment's split share. Taken: per-level effects weighted by treated rows (8 to
   2,561 rows per level, verified); the treatment column's split share reported. ACIC is stated as a
   guard on one covariate set.
8. Pointwise bias. Taken: test points fixed across replications, bias as the root mean square over points
   of each point's mean error (Morris, White and Crowther 2019).
9. Citations and Strobl's design. Taken: Breiman et al. 1984, Loh and Shih 1997, Kim and Loh 2001,
   Hothorn, Hornik and Zeileis 2006, Bleich et al. 2014, Benavoli et al. 2016 and 2017; Strobl's design
   verified (M(2), M(4), M(10), M(20), n = 120, relevance 0.05 to 0.2) and its bias named as a search
   mechanism.
10. Order and scope. Taken: real blocks first; ACIC and the long chains only on an R1 pass; flexBART at
    n = 1000 only; Friedman at n = 500 only; the Strobl null cut to 20 replications since p5 settled it;
    R on the real regression block and ACIC only.

Also changed: the equal-time rule no longer has a 10 percent trigger; any passing candidate slower than A
is re-read at matched time.

## Amendments after the re-check

The critique's re-check of the first amendment (verdict: run after amendments) made ten more; each was
checked first, by rerunning its probe (scratch/fsp-crit/probes/bsr-check.R, which reproduced every number
it quotes) and a new one (probe bsr2).

1. A per-dataset harm guard. Done: R1 now fails to a reported trade-off when any dataset's correlated
   Bayesian t-test (Benavoli et al. 2017) puts more than 0.95 on a practical loss. Checked: one dataset
   losing 49 percent beside 13 gaining 2 percent, which the signed-rank test alone passed 0.62 of the
   time, now never passes; the guard costs little power (a uniform 2 percent gain, 0.88 to 0.87 precise).
2. ACIC's settings fixed or sampled. Done: fixed, the competition's whole design; the standard error and
   the pilot now use within-setting variance alone, so the count does shrink with replications. Taken as
   sampled, a true effect varying across settings by more than about 0.035 would put the bound out of
   reach at any count (the critic's figure), and the estimand would be a population of settings the
   competition never defined.
3. Coverage margins and the side of the bound. Done: coverage measures take 0.02, the maintainer's bar for
   a coverage gap (dec-B455), with that study's noise correction; the bound is one-sided 97.5 percent
   (ICH E9). The 1-point margin carried from error had no argument and was out of the budget's reach.
4. A bias guard for SATT, or drop "bias". Done the second way: SATT's bias is inside the RMSE guard and is
   reported; the summary no longer says the guard covers bias.
5. R1's power and false-gain rate at both error scenarios. Done in the sizing table; the multiplicity
   statement now rests on the noisy rate (1.5 percent regression, 0.7 percent binary per candidate; at
   most about 7 percent family-wise by the union bound), not on "0 of 400".
6. "Equivalent" as outcome error. Done in the rule and the summary.
7. The span floor on the binary block. Done (adult's 41-level factor and cleveland, infert and birthwt,
   all under 6 levels, are expected to meet it; step 1 checks).
8. The indicator arm only through the pilot. Done: R runs on ACIC replications 1 to 10 only; the pass
   path drops about 10 core-hours.
9. M for any level frequencies, derived here. Done; the cancellation is written out under Candidates.
10. The open call's background. Done: it names the prior mechanism, the permutation null (Bleich et al.
    2014) and that option (a) saves about 9 core-hours.

Numbers that changed: the pass path is 20 to 41 core-hours (25 to 46 with DART), from 20 to 51 (25 to
56); the critic's re-check expected about 50 with the indicator arm run throughout. R1's quoted power at
a uniform 2 percent gain is 0.87 at precise errors and 0.39 at noisy ones, where it was 0.86 alone.

## Not verified

- The x86 box's speed against the laptop (the 1.5 factor); the smoke re-prices.
- That all fourteen regression datasets can be fetched and matched to Deshpande's Table 1; his
  preprocessing is published only as the table's dimensions. That fuelEcon's near-zero SMSE has the same
  cause as Caschool's (a response defined from predictors) is not checked; it is left out on the published
  scores.
- That mushroom is separated almost exactly by single predictors is recalled (odor), not run here.
- The per-fit spread of ACIC's effect measures and the arms' coverage discordance, which set R2's
  replication count; the pilot measures them. That the 77 settings were built as a grid rather than drawn
  is read from the competition's design, not rechecked here.
- Cited from knowledge, not read here: Breiman et al. 1984, Loh and Shih 1997, Kim and Loh 2001, Hothorn et
  al. 2006, Bleich et al. 2014, Benavoli et al. 2016, Morris et al. 2019, Berger 1982, Wittes and Brittain
  1990, ICH E9 (1998), Nadeau and Bengio 2003, Dorie et al. 2019's scoring, Linero 2018. Read: CGM 2010,
  Deshpande 2024, Linero and Yang 2018, Strobl et al. 2007, Benavoli et al. 2017, and the package sources in
  the table.
- Time on factors of more than 63 levels in a 20,000-row dataset beyond the critic's 84-level probe.
