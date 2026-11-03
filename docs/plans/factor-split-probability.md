# factor-split-probability: should a factor's prior split probability grow with its number of levels?

Status: PLANNED 2026-11-03, before its blind critique. Research item of TODO "factor encoding" (dec-B78),
scheduled between the merge to main and the 1.0-0 submission (VD 2026-09-28). Nothing has run beyond
the sizing probes below.

agent: the blind critique (opus) follows this draft. One opus implementer writes the harness, runs the
smoke on the laptop and the main run on the bench box, and drafts the report; one opus review of the
harness and smoke before the main run, since the harness carries the statistics; two independent
referees (opus) on the report before it reaches the maintainer.
rng: NEUTRAL. No package code changes: every candidate rule is a `split.probs` vector the shipped
package already takes, and the BayesTree reference arm is `factors = "indicators"`. Only the
conditional DART stage needs a scratch build (about 30 lines, never landed). A change of default the
study supports is a separate, draw-changing slice with its own plan.
window: after the merge to main, before 1.0-0 (dec-B78, VD 2026-09-28).
budget: harness about 900 lines of R (benchmarks/R/factor-split-probability.R) reusing the house's
binary dataset loaders; compute about 26 box core-hours expected, 35 with the effects confirmation
block, stop at 55 ([Sizing](#sizing)); report about 5 pages of text plus figures.

## Summary for the maintainer

**The question.** When a tree picks a variable to split on, should a factor with many levels be picked
more often than a numeric column, in proportion to its number of levels or some function of it?

**Why it is open.** The BayesTree-style coding turns a factor with K levels into K indicator columns,
each picked as often as a numeric column, so a 16-level factor was picked 16 times as often as a
numeric column. That was never argued for; it came with the coding. Since 1.0-0 a factor is one column
split on any subset of its levels (dec-B78), so it is picked as often as any other column. You asked
whether the old weighting should come back on the new coding.

**What others do.** Every package that sets a factor's weight on purpose gives it one column's weight:
the original Bayesian CART prior (pick a variable uniformly), bartMachine (which divides a factor's
weight among its indicators, since its version 1.0.2), stochtree (the same), SoftBart's formula
interface and Linero's grouping prior (one share per factor, split among its indicators), flexBART
(one column, subset splits, as dbarts now) and pymc-bart (one column). Only BayesTree and the BART
package weight a factor by its level count, and both only because they code factors as indicators;
the BART package computes a grouping of a factor's indicators and its sampler does not read it.

**What the probes found.** Two effects pull opposite ways, and both are large enough to see in ten
replications.

- Where a factor drives the function (the four published flexBART designs, one 10-level factor and ten
  numeric columns), weighting by level count cut the error of the fitted function by 0 to 8 percent.
- Where factors carry no signal (a standard function plus factors of 2 to 50 levels that do nothing),
  the same weighting raised the error by about 50 percent and widened intervals, because the trees
  spent nearly a third of their splits on a useless 50-level factor.
- On three small real datasets the differences were within noise.
- With no signal, the share of splits each predictor takes after fitting equals the share the prior
  gave it. So a weighting rule biases split-count importance toward many-level factors by exactly its
  weight, which is the defect the random-forest literature names in importance measures.

So which rule is better depends on which kind of data users meet, and the study reads that on real
data.

**Candidates.** Today's rule (one column's weight); the level count (BayesTree's implicit rule: K for
K of 3 or more, 1 for a binary factor); the level count less one (the factor's degrees of freedom);
the square root of that (the midpoint of that family, which no source argues, included as the BART
paper weighed its defaults against a grid of settings); and, as a reference only, the BayesTree
indicator coding.
Weighting by the log of the number of possible subset splits is the level count less one to within 0.2
for 4 or more levels and is not fitted separately.

**Data, mostly chosen by others for the purpose.**
- Sixteen real regression datasets with factors, the benchmark the flexBART paper assembled for
  categorical predictors in BART (2 to 84 levels), scored on held-out rows as the BART paper scored its
  defaults; and the house's real binary datasets that have a factor of 3 or more levels (credit,
  census income, heart disease and others), scored by held-out log score.
- The 2016 causal inference competition's data: real covariates with factors of 6, 16 and 5 levels and
  published surfaces with known treatment effects, scored as the competition scored them (average
  effect on the treated: bias, error, interval coverage; individual effects).
- Published simulated designs that explain the mechanism: the four flexBART designs, and the
  no-signal and one-signal designs of Strobl et al. (2007), the standard test of whether an importance
  measure prefers many-level factors; plus the BART paper's test function with null factors.

**What decides, fixed before any fit.** Today's rule stays unless a candidate is better on the real
datasets (median across datasets better, and better by a signed-rank test across datasets, the
standard way to compare methods over many datasets) and loses nothing on the competition's treatment
effects (error, bias and interval calibration, over- and under-coverage both counting). The designed
cells explain the result and do not decide it. If a weighted rule wins, a DART check runs before any
recommendation.

**Your call** (below): whether the bias in split-count importance a weighted rule brings by
construction is a reason on its own to keep today's rule. Recommended: no; report it, and say it in
the manual if a weighted rule wins, since the estimates are what users act on.

**Cost.** About 26 core-hours on the box, 35 if the effects need their confirmation block; stop at 55.
The result comes to you as a short PDF report after two independent referees.

## Goal

A research report (a PDF of at most five pages of text plus figures, self contained) that answers
whether the default `split.probs` for a single-column factor should depend on its level count, argued
as the BART literature argues for its tree-prior defaults, refereed twice, with a recommendation the
maintainer can rule on. Nothing in the package changes under this plan.

## Context

- dec-B78: "Weighting a factor column's split mass by its number of levels is a research item, moved on
  2026-09-28 from after the release to between the merge and the submission ("I would actually move up
  the research on split probability to post-merge, pre-release"). The maintainer chose "Subset splits
  as the modern default."" The encoding is ruled; only the weighting is studied.
- dec-B559: the study "hold[s] up to academic scrutiny", comes as "small research reports (pdf, at most
  5 pages of text with as many graphs as needed, self contained and easily read)", and imitates "how
  other people argue for default priors, not inventing things ourselves." docs/plans/README.md,
  Studies, adds the dec-B558 questions: datasets people meet, natural estimands, RMSE, bias and
  heterogeneous effects.
- dec-B512 (revising dec-B116, dec-B246, dec-B252): "split.probs takes a vector of split probabilities
  or dart(), which is DART, a prior on which variable a split uses". A changed default is the value
  `split.probs` takes when unstated, at each door that takes it.
- dec-B338: "split.probs weights how often each predictor is picked for a split"; negative, non-finite
  and all-zero vectors are refused. dec-B547: a named vector reads by column or term name with an
  optional `.default`.
- dec-A09, dec-A115: an ordered factor splits between its level codes, as a numeric column does, from
  its declared levels. Out of scope: it is weighted as a numeric column is.
- dec-B370, dec-A26: wide factors are stored sparse and split by subset like any factor; the rule
  applies to them unchanged.
- The engine: `CGMTreePrior::drawSplitVariable` picks a column among those available at a node with
  probability proportional to its `split.probs` entry (uniform when none), and
  `CGMTreePrior::drawRuleForVariable` then draws a factor's rule uniformly over the left-right patterns
  of its reachable levels that leave neither side empty. `DartPrior` draws the probabilities from a
  symmetric Dirichlet, alpha / p per column plus the split counts.
- The help (man/bart.Rd, `factors`) already says subset splits "give a factor the prior split
  probability of one predictor rather than one per level".

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

The recursive-partitioning literature on variable selection bias (here Strobl, Boulesteix, Zeileis
and Hothorn 2007, BMC Bioinformatics) holds that "in the null case ... a sensible variable importance
measure should not prefer any one predictor variable over any other", and tests it with one continuous
predictor and four factors of 2 to 20 categories at n = 120. Deshpande (2024) shows one-hot coding
cannot pool levels and that uniform subset splits are "biased towards balanced partitions and ... generally
unable to isolate a single outlying level" (his DGP2); his benchmark compares flexBART's coding and
weight with BART's coding and weight together, so it cannot say which of the two did the work. This
study separates them.

## Probes (laptop, shipped build, 2 cores, 2 chains; scratch/fsp-plan/probes)

| probe | design | finding |
|---|---|---|
| p1 | flexBART DGP 1-4, n = 1000, 10 reps | log RMSE ratio, level count vs uniform: -0.07, -0.04, -0.08, 0.00 (sd per rep 0.06 to 0.07); indicators vs uniform +0.11, -0.09, +0.11, +0.06 (Deshpande's pattern); coverage of f, paired difference -0.01 (sd 0.02); posterior split share of the factor 0.24 uniform, 0.51 weighted, against prior 0.09 and 0.50 |
| p2 | CGM's Friedman function on 5 of 10 numerics plus null factors of 2, 4, 10, 20, 50 levels, n = 500, sigma 1, 10 reps | log RMSE ratio +0.41 (sd 0.13) for level count, +0.36 level count less one, +0.33 indicators; coverage 0.934 uniform, 0.898 level count; the null 50-level factor took 0.30 of the splits under the level count, 0.026 under uniform |
| p3 | real: ggplot2 mpg (factors of 3 to 38 levels), carData Salaries, MASS Cars93; 20 splits of 75/25 | log RMSE ratio level count vs uniform +0.04, -0.003, +0.05, per-split sd 0.20, 0.007, 0.16; predictive 90 percent coverage within 0.03 across arms |
| p4 | ACIC 2016 covariates, a toy surface, default `bart` (4 chains, 2 threads) | 5.6 s per fit (about 11 core-seconds); the level count arm and uniform within 5 percent in time, indicators 30 percent faster |
| p5 | Strobl null, gaussian response, n = 120, 20 reps | posterior split share equals the prior share: uniform 0.20 each; level count 0.03, 0.03, 0.11, 0.28, 0.55 against prior 0.03, 0.03, 0.11, 0.28, 0.56; sd per rep under 0.008 |

## Candidates

Weights apply to an unordered factor's column when `split.probs` is unstated; numeric, logical and
ordered-factor columns keep weight 1; K is the number of levels with training rows.

- A, uniform (today): w = 1. A factor is one predictor. Split counts compare fairly across columns
  under no signal (p5). A many-level factor that matters is proposed 1/p of the time.
- B, indicator count: w = K for K >= 3, 1 for K = 2. BayesTree's and 0.9-34's prior mass on the new
  coding. A 50-level factor beside 10 numeric columns takes 83 percent of the prior splits, and split
  counts favor it by that much with no signal.
- C, degrees of freedom: w = K - 1. B less one; the two differ materially only at K = 3 or 4. No source
  argues it as a split prior; it is the factor's parameter count in a linear model and agrees with
  every rule at K = 2.
- S, the family's midpoint: w = sqrt(K - 1). Not argued anywhere; included because CGM 2010 (sec 5.1)
  weighed its defaults against a grid of settings, and because B and C share a direction (the whole
  level count) the probes show can hurt; S tests whether part of the weight keeps the gain.
- R, reference only: `factors = "indicators"`, the BayesTree coding and weight. Not a candidate (dec-B78
  rules the coding); it lets the report say what 0.9-34 users had.
- Considered and not fitted: the log of the number of distinct subset splits, log2(2^(K-1) - 1), which
  is within 0.2 of C for K >= 4 (1.58 against 2 at K = 3, 2.81 against 3 at K = 4); a node-local weight
  from the levels reaching the node (what indicator coding does at depth, since a separated level's
  indicator has no split left), an engine change reached only if B or C wins; and BART's K = 2 weight
  of 2, an artifact of its coding that its documentation contradicts ("one dummy is created if q=2").

**With a stated `split.probs`.** A stated vector is literal (dec-B338) and no rule applies. How a
named vector's `.default` reads under a weighted default is an implementation-time call, not this
study's.

**With DART.** `dart()` takes no vector; its Dirichlet is symmetric, alpha / p per column. A weighted
default carries to DART as its prior mean, alpha * w_j / sum(w) per column, a non-uniform base measure
for its Dirichlet. Since the split counts dominate the Dirichlet once trees use a column, the weight is
expected to matter less there. This is read only if a weighted candidate passes the decision
rule (stage 6), on a scratch build of that base measure.

**With forest column restrictions and interaction constraints.** The weights renormalize over the
columns a forest or node may use, as today.

**Mixing against prior.** A weight changes the prior and how often the sampler proposes the factor at
once. At equal running time the study measures both, as users get them; the long-chain cell (below)
says which one a gain comes from. A gain that is mixing alone would point to a proposal change that
keeps the uniform prior, a separate design, not to a new default.

## Estimands and scores

Every summary is computed per draw and averaged last (posterior means of the quantity asked for).
Intervals are central 95 percent posterior intervals for f and effects, 90 percent predictive intervals
for new outcomes, scored as distance from nominal, over- and under-coverage both counting; width is
reported beside coverage.

- f(x) at held-out x where the truth is known: RMSE, mean signed error (bias), coverage and width.
- Level means, the average of f over the rows at each level, the quantity an analyst reports for a
  factor: RMSE, bias and coverage over levels.
- Treatment effects (ACIC 2016), as Dorie, Hill, Shalit, Scott and Cervone (2019, Statistical Science)
  scored the competition: the average effect on the treated (bias, RMSE, coverage and width of its
  interval), individual effects (PEHE, the RMSE of individual effects, and coverage of their
  intervals), and effects averaged within each level of the 16-level factor, the heterogeneous effect
  over groups.
- Real outcomes: held-out RMSE (CGM 2010 sec 5.1; Deshpande's SMSE is the same ratio up to a constant),
  summarized as CGM's relative RMSE (RMSE over the best arm's on that split); coverage of predictive
  intervals; for binary outcomes, held-out log score and Brier score, as the house's binary studies
  scored them.
- Variable importance, as the literature scores it: the share of splits each predictor takes (CGM 2010
  sec 5.2 reads usage frequency on the Friedman data) and the share of draws using it at all; under
  no signal, against equal shares (Strobl et al. 2007); under one signal, the informative predictor's
  rank and share.
- Time per fit, for the equal-time read.

**Equal time.** Fits run at the shipped defaults (4 chains, 500 burn-in, 500 draws, 75 trees). If a
candidate's median time per fit in a component exceeds uniform's by more than 10 percent, uniform is
rerun in that component with burn-in and draws scaled by the time ratio and compared at that length.

## Data

Every source but two was published by others for the purpose it serves here: the house's binary set
predates this question, and the Friedman cell with null factors joins two published designs. None was
picked after seeing a fit.

1. **Real regression datasets with factors**: the sixteen of Deshpande (2024), Table 1, the benchmark
   assembled for categorical predictors in BART (abalone, ais, Alcohol, amenity, attend, cane,
   Caschool, cpu, Engel, fuelEcon, Insur, Medicare, mpg, servo, spouse, strike; n from 167 to 23,971;
   1 to 8 factors of 2 to 84 levels), with his sources (UCI, the Journal of Applied Econometrics data
   archive, R packages, fueleconomy.gov, StatLib, OzDASL). Step 1 fetches each, caches it with a
   recorded sha256 as benchmarks/R/uci-binary.R does, and checks its dimensions and level counts
   against Table 1; a dataset that cannot be fetched or matched is recorded and left out, with no
   substitute. 20 random 75/25 splits each (CGM 2010 used 20; Deshpande 50; the sizing below shows
   more splits buy little), test rows at a level absent from training dropped and counted (dbarts
   refuses an unseen level, dec-A85), the same rows for every arm.
2. **Real binary datasets with factors**: those of the house's 22 (benchmarks/R/binary-datasets.R,
   uci-binary.R; assembled for the binary k studies before this question) with an unordered factor of
   3 or more levels: adult, bank, mushroom, german, tictactoe, creditapproval, cleveland, infert,
   birthwt (step 1 confirms the list from the loaders' columns). The house's split and cap: 80/20, at
   most 4,000 training rows past 5,000 rows. 20 splits each.
3. **ACIC 2016** (Dorie et al. 2019; aciccomp2016): 4,802 rows of real covariates from a perinatal
   study, 58 columns with factors of 6, 16 and 5 levels; all 77 published settings, replications 1 to
   5, the competition's known potential outcomes as truth. A fit puts the treatment as a predictor and
   reads effects from the two counterfactual predictions (Hill 2011, JCGS).
4. **Designed, to explain the mechanism (do not decide)**:
   - flexBART DGP 1-4 (Deshpande 2024, sec 4.1): ten uniform numerics and one 10-level factor, n = 1000
     and 5000, 25 replications, 500 test points; the factor drives the function.
   - Strobl et al. (2007): X1 continuous, X2 to X5 factors of 2, 4, 10 and 20 categories, n = 120; the
     null case and the power case (the response depends on the binary X2, at three strengths), binary
     response as theirs and a gaussian null; 50 replications.
   - CGM 2010's Friedman function (their sec 5.2) on 5 of 10 uniform numerics, sigma 1, plus four null
     factors at Strobl's level counts and one of 50 levels (site- or state-like), n = 500 and 2000, 25
     replications; the factors carry nothing.
   - Long chains: flexBART DGP 1-4 at n = 1000, arms A and B, 10 replications at four times the
     default length, to separate a mixing gain from a prior gain.

## The decision rule (fixed before any fit)

Sources: CGM 2010 sec 5.1 (defaults argued by out-of-sample RMSE over a benchmark of real datasets,
read as relative RMSE across datasets and splits); Deshpande (2024) sec 5.1 (per-dataset paired tests
with a Bonferroni correction, and the count of datasets won); Demsar (2006, JMLR) for comparing methods
over many datasets by a signed-rank test across datasets; Dorie et al. (2019) for the causal scores.

For candidate c and real dataset j, L_cj is the mean over splits of log(RMSE_c / RMSE_A) (for binary
data, of the held-out log loss ratio).

- **R1, real outcomes.** In the regression block and in the binary block separately: the median over
  datasets of L_cj is below 0 in both, and in at least one block the Wilcoxon signed-rank test over
  datasets is significant at 5 percent two-sided after Holm's correction over the three candidates,
  with neither block significant against c. Per-dataset paired tests (Nadeau and Bengio's correction
  for overlapping training sets, Bonferroni over datasets and candidates) and CGM's relative-RMSE
  boxplots are reported, not decisive.
- **R2, treatment effects (ACIC 2016).** Paired over settings and replications, c is not significantly
  worse (5 percent two-sided, Bonferroni over the measures) than A in RMSE of the average effect on the
  treated, the distance of its coverage from 0.95, PEHE, the distance from 0.95 of the coverage of
  individual effects and of f, or the RMSE of the per-level effects. R2 is read for a candidate that
  passed R1 only after replications 6 to 10 have run for it (the confirmation block).
- **Verdict.** A candidate passing R1 and R2 is recommended; among several, the one with the lowest
  median L. Otherwise today's rule stays. A candidate passing R1 and failing R2 (or the reverse) is
  reported as a trade-off and today's rule stays.
- The designed cells, variable importance and the long chains are reported beside the verdict to
  explain it. Whether the importance bias of a weighted rule is decisive on its own is the open call.

## Sizing

Smallest decision-changing difference: R1 turns on signs and ranks across datasets, so per-dataset
precision matters only enough to keep a dataset's sign from being noise; R2 on a significant loss.

| component | question | size | box core-h | precision reached (from probes) | minimum that resolves | sequential stop |
|---|---|---|---|---|---|---|
| real regression | R1 | 16 x 20 splits x 5 arms = 1,600 fits | 8.7 | per-dataset SE of L under the overlap correction: 0.62 x per-split sd (0.004 on Salaries-like data, 0.12 on mpg-like); at 50 splits 0.59 x, so more splits barely help | 20 splits (CGM's count) | none: extra splits do not move the corrected SE |
| real binary | R1 | 9 x 20 x 5 = 900 | 3.0 | not probed; the house's binary studies' split variance | 20 | none |
| ACIC 2016 | R2 | 77 x 5 x 5 arms = 1,925 | 8.8 | coverage SE about 0.011 per arm over 385 fits, a paired difference about 0.006; PEHE log ratio SE about 0.01 at a per-fit sd of 0.2 | 5 reps for the read, 10 to confirm | replications 6 to 10 only for a candidate passing R1 (+8.8) |
| flexBART DGPs | explain | 4 x 2 n x 25 x 5 = 1,000 | 3.3 | log RMSE ratio SE 0.014 (p1 sd 0.07); paired coverage SE 0.004; level means 0.044 | 25 resolves 3 percent and 0.01 coverage | stop at 25 |
| Strobl | importance | 4 cases x 2 families x 50 x 5 = 2,000 | 1.2 | split share SE about 0.001 (p5 sd 0.006) | 20 would do; 50 for the power case's ranks | stop at 50 |
| Friedman + null factors | explain | 2 n x 25 x 5 = 250 | 0.5 | log RMSE ratio SE 0.03 (p2 sd 0.13) on an effect of 0.4 | 10 | stop at 25 |
| long chains | mixing or prior | 4 x 10 x 2 arms, 4x length | 0.4 | the p1 gap of 0.07 at SE 0.02 | 10 | none |
| DART (conditional) | does a weighted base measure hold | the passing candidate and A, uniform-DART and weighted-DART, real and ACIC blocks at half size | about 5 | as above at half size | | only if a candidate passes |

Box core-hours are laptop core-seconds times 1.5, an unmeasured ratio the smoke re-prices; fits at the
default 4 chains. Expected 26; 35 with the ACIC confirmation; about 40 with the DART stage; stop at 55.

## Steps

1. Branch work only. Harness benchmarks/R/factor-split-probability.R: loaders for the sixteen
   regression datasets (fetch, sha256, dimension check against Deshpande's Table 1), the binary list
   from the house's loaders, ACIC 2016 through aciccomp2016, the designed generators, the arms as
   `split.probs` vectors computed from each training set's level counts, one block per invocation
   writing one rds, and a summarize mode that computes R1, R2 and the reported tables. Verify: every
   loader's dimensions printed against the table; the prior split share of every factor under every
   arm for every dataset (no fits), which goes into the report as its first figure.
2. Smoke on the laptop at 2 cores: one split or replication per block per arm; time per fit per
   component, re-pricing the table above; the equal-time check.
3. Opus review of the harness and smoke (the statistics: pairing, per-draw summaries, the overlap
   correction, the rank test, the truth extraction from ACIC). Fix round.
4. Pre-registration: commit the plan, any amendment from step 3, and the harness at one hash before
   the main run; amendments after that are recorded as such.
5. Main run on the box: real blocks, ACIC replications 1 to 5, designed cells, Strobl, long chains.
   Read R1. For a candidate passing R1, ACIC replications 6 to 10, then R2.
6. Conditional: if a candidate passes R1 and R2, the DART stage on a scratch build whose base measure is
   alpha * w_j / sum(w), and the node-local variant considered above is named as follow-up work.
7. Report: R Markdown source docs/design/factor-split-probability-report.Rmd rendering every figure
   from the stored outputs; at most five pages of text. Outline: (1) the question and why it is open;
   (2) what the literature and packages do (the table above); (3) the candidates and what each means
   for a user, with the prior share figure; (4) data and estimands, how each was chosen; (5) results:
   CGM-style relative RMSE boxplots by dataset, the binary block, ACIC's effect scores, designed cells
   (error and coverage by arm), the Strobl null and power shares (their Figure 1's form), the long
   chains; (6) recommendation and limits. A finding note docs/design/factor-split-probability.md
   records the verdict with a link to the report.
8. Referees: two independent opus referees, each on a neutral brief (scratch/fsp-peer-common.md,
   modeled on scratch/nbs-peer-common.md: referee it as a methods-journal submission, owe it no
   deference, verify before asserting, read-only access to the stored outputs, small simulations of
   their own allowed, verdict first, ranked findings, the smallest addition that would change their
   verdict priced in core-hours). The brief names no preferred outcome. The report is revised to
   answer them; the maintainer receives the report, both referee reports and the response.
9. Records: plan Status and a Landing note, the finding, TODO "factor encoding", the INDEX rows; the
   ruling goes to the ledger when the maintainer makes it.

## Constraints

- No package code changes; NEUTRAL. A recommended change is its own plan and slice (draws change).
- At most 2 cores on the laptop; the box per its own limits; box folders deleted at landing, study
  data archived to the laptop first.
- No reference chains as truth; real outcomes scored on held-out rows only.
- The decision rule above is frozen at step 4; any amendment after a fit is scored is labeled as such.
- Out of scope: the coding (dec-B78), ordered factors (dec-A09), `bartBT`, numeric columns' weights,
  the subset-split rule's own distribution over partitions (Deshpande's DGP2 point), a proposal change
  that keeps the prior.

## Verification

- Step 1: `Rscript benchmarks/R/factor-split-probability.R datasets` prints each dataset's rows and
  factor level counts beside Deshpande's Table 1, with every mismatch marked; `prior` prints the
  prior share table.
- Step 2: `Rscript benchmarks/R/factor-split-probability.R smoke DIR` completes every block, and
  `summarize DIR` refuses to mix smoke and main outputs.
- Step 5: `summarize DIR` prints R1 and R2 with the counts of fits per cell, no cell short.
- Lint and air on the harness; `tools/check-doc-freshness.R`.

## Open calls for the maintainer

1. **Does split-count importance decide on its own?** Background: with no signal, the share of splits
   each predictor takes after fitting equals the prior share (probe p5), so any weighted rule makes
   split counts favor many-level factors by exactly its weight; the random-forest literature calls that
   a defect in an importance measure. dbarts reports split counts, and users read them as importance.
   Options: (a) importance decides alongside the estimates: every weighted rule fails it by
   construction, so today's rule stands without the study, saving its cost; (b) the estimates decide
   (R1 and R2), importance is reported, and a winning weighted rule ships with the manual saying split
   counts then favor many-level factors in proportion to their weights; no extra cost; (c) the
   estimates decide, and a winning rule must also not lose the informative predictor in the Strobl
   power case; no extra cost, a rule that may fail on a 120-row design. Recommended (b): estimates and
   intervals are what users act on (your stated priority, 2026-10-17, that f and treatment effects come
   first), and a user who wants
   even counts can state `split.probs`.

## Not verified

- The x86 box's speed against the laptop (the 1.5 factor); the smoke re-prices.
- That every one of Deshpande's sixteen datasets can be fetched and matched to his Table 1; his
  preprocessing is not published beyond the table's dimensions.
- The level counts in Strobl et al.'s Table 1 (the text, read, gives 2 to 20 categories and n = 120;
  the 2, 4, 10, 20 split is recalled, not read).
- Nadeau and Bengio (2003) and Demsar (2006) are cited from knowledge, not read here.
- Which of the house's binary datasets carry an unordered factor of 3 or more levels; step 1 confirms.
- Time on wide factors (more than 63 levels, the pooled masks) and on the 20,000-row datasets.
- Linero (2018), the DART paper, was not read here; its Dirichlet is read from dbarts' and BART's code.
