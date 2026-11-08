# study-firming: the first tier of landed studies brought to published practice

Status: planned (2026-11-07), revised after referee round 1. Nothing has run. Five studies, each separable, in the
order the maintainer set (dec-B637): the swap move, the logistic k default, the negative binomial default,
stan4bart's initial fit with the starting sigma, and the monotone defaults. Each study's design and readings are
written here before any fit it reads. No reading carries a decision rule: each report gives its readings with
their intervals and recommends, as dec-B625 did, and the maintainer decides from a report of at most 5 pages of
text (dec-B559).

agent: per study, one opus implementer for the harness and the runs on the bench box, one review of the harness
and its smoke (the gate-only tier of [Process by risk](README.md#process-by-risk)), and two independent referees
on the report. The rereads need no implementer: one agent on the laptop.
rng: the studies change no package code. Every run is on one 1.0 library in shipped mode; an arm is a user-level
setting, never a build.
window: before 1.0-0. The logistic k study's Part D waits for k-move-pricing's ruling on the move (dec-B634).
Studies run one at a time on the box.
budget: about 31 core-hours expected over the five, about 46.5 at most, before the options in
[For the maintainer](#for-the-maintainer) ([Sizing](#sizing)); each study's harness a few hundred lines of R,
much of it reused from the landed studies'.

## Goal

After this plan runs, each of five landed defaults rests on a study that argues as published work of its kind
argues: its question in a user's terms, its estimands, comparators, datasets and run lengths taken from a named
published comparison, its readings shown as those reports show them, and every choice without a published model
named as the project's own (dec-B633). Each study ends in a recommendation the maintainer rules on: keep the
default, change it, or reword what the documentation claims for it.

## Context

- The standard: [Studies](README.md#studies) (dec-B559, dec-B626, dec-B633, dec-B637). Readings without a
  published bar are reported without a bar, with their intervals, as dec-B625 did.
- The landed studies: [swap-removal.md](../design/swap-removal.md) with sections 13 and 14 of
  [tree-mixing-proposals.md](../design/tree-mixing-proposals.md);
  [binary-k-fixed-or-drawn.md](../design/binary-k-fixed-or-drawn.md);
  [nbinom-default-study.md](../design/nbinom-default-study.md);
  [stan4bart-response-range.md](../design/stan4bart-response-range.md);
  [starting-sigma-sensitivity.md](../design/starting-sigma-sensitivity.md), its section on the forests' bases;
  [monotone-exact-birth-death.md](monotone-exact-birth-death.md), its feel and tree-prior studies.
- What 0.9-34 did, run on the released package: the tree-move mixture was birth/death 0.5, swap 0.1, change 0.4
  (birth 0.5 within birth/death); a binary response drew k under the improper chi(1.25, Inf); there was no
  logistic or negative binomial family, no monotone constraint and no basis. stan4bart 0.0-13, the release that
  paired with 0.9-34, no longer loads beside today's libraries (its compiled threading library does not match);
  its source fits the initial mixed model without the forest's columns and hands its fitted values to the forest
  as an offset.
- k-move-pricing (on its branch) owns the logistic sampler question and a held-out list of datasets. Nothing here
  reads, names or fits a dataset from that list (dec-B635).

## Common to every study

- Defaults as a user gets them: `bart()` with 75 trees, 4 chains, 500 burn-in and 500 kept per chain, unless a
  reading copies a published run length, which is then named.
- Paired arms: every arm of a replication shares its data seed and its sampler seed, as Morris, White and
  Crowther (2019, section 5.4) advise, which lowers the Monte Carlo SE of a comparison between methods.
- Every reading is a mean over replications with its Monte Carlo SE and 95 percent interval (Morris et al.,
  Table 6 and section 6.2); a difference between arms is a paired difference with its interval. No rule reads an
  interval, so none is adjusted for the number shown; each report says how many it shows, of which about one in
  twenty excludes zero by chance at 95 percent (the project's own; alternative: Bonferroni intervals, or
  Holm-adjusted p-values beside the plain intervals).
- Reporting: each report states its design in Morris et al.'s ADEMP order (aims, data-generating mechanisms,
  estimands, methods, performance measures; their Table 1 and section 6.1) and recommends from the intervals,
  reading with no bar as the published report it follows reads its own: Tan et al. (2024, L.6), "does not
  substantially affect"; CGM (2010, 5.1), "arguably second best"; mBART (5.3), "good reliability".
- Sizing (Morris et al., section 5.3): for coverage of one interval per replication (an average effect, a fixed
  effect), eq. (1), n = p(1 - p) / MCSE^2 at p the nominal level: at 95 percent the interval's half-width is 0.085
  at 25 replications, 0.060 at 50, 0.043 at 100, 0.030 at 200 and 0.020 at 456. For coverage averaged over test
  points, a replication's coverage is a mean of indicators, whose variance is at most p(1 - p), so eq. (1) bounds
  the size; the size used is Table 6's form for a mean, MCSE = s / sqrt(n), with s the per-replication spread from
  the stored runs, as Morris et al. size from an initial run (section 7.1), rechecked at each smoke. For a paired
  difference in single-interval coverage the variance is d minus the squared difference, d the share of
  replications where exactly one arm covers: at 200 replications the half-width is 0.031 for d = 0.05 and 0.044
  for d = 0.10. The targets each size aims at are the project's own, set at the landed studies' resolution;
  Morris et al. give the method, not the targets.
- Machine: the x86 bench box for every run, the laptop for every reread. Costs below are from single-thread
  timings on the laptop under a load of 8 to 11 and are rechecked in each harness's smoke.
- Stop rule (the project's own, as pricing's is): two looks, at half and full size. A component stops at the
  first look once every interval it reports is inside its target in the sizing table. A cell still wider than its
  target at full size is reported at the width reached, with the replications eq. (1) or Table 6 says would meet
  it. No component runs past 1.5 times its expected cost: it stops and reports what it has.
- Each report: at most 5 pages of text; each table or figure tied to the published table or figure it copies;
  the study's own choices listed with their alternatives; a section on what it does not settle.

## A5. The swap move dropped from the default mixture

**Question.** A user who never sets `proposal.probs` gets a mixture without the swap move: birth/death 0.6,
change 0.4. Does that change the accuracy, the interval coverage or the mixing of what they estimate, against the
0.9-34 mixture? The news entry says "dropping it costs nothing measurable in accuracy, coverage or mixing".

**Published reports.**
- Tan, Ronen, Saarinen and Yu (2024), appendix L.6, "Experiment 7: Restricted vs. full move set": the nearest
  published contrast, and a larger one: the full mixture (grow .25, prune .25, change .4, swap .1) against grow
  and prune alone ("constraining the proposals to just grow and prune"), which drops change as well as swap. 200
  trees, 8 chains, 1000 burn-in and 10,000 kept, 25 replicates; training sizes 200, 500, 1000 and 10,000 drawn
  from four PMLB datasets (Breast Tumor, Echo Months, Satellite Image, California Housing; a fixed 10 percent test
  set) and two synthetic designs (Low-Dimensional Smooth, Piecewise Linear). Readings (L.1): on real data RMSE
  against y and coverage of 95 percent posterior predictive intervals for y (eqs. 127 and 128); on synthetic data
  RMSE and coverage of credible intervals for f (eq. 129); Gelman-Rubin R-hat in its classical form across
  chains; bars of 1.96 standard errors. Verdict, read with no bar: "restricting the move set does not
  substantially affect R-hat, coverage, or RMSE."
- Chipman, George and McCulloch (1998), single tree: "omission of the SWAP step inhibits mixing".
- Pratola (2016): effective sample size of sigma and per second (Table 2), run times (Table 3).
- He and Hahn (2023), Table 4: coverage, interval length, RMSE and time over 100 replications.
- Vehtari et al. (2021): 4 chains, rank-normalized split R-hat below 1.01 and bulk ESS of at least 400.

**Design.** Two arms: the shipped mixture and the 0.9-34 mixture, both set through `dbartsControl` on one library.
- C1: the landed study's five designs (Friedman at two noise levels, a depth-3 interaction, 100 predictors, a
  probit), 25 replications at the defaults.
- C2: Tan's six datasets at n 500 and 1000, two of Tan's sizes, 25 replications at the defaults, Tan's test sets.
- C3: Tan's settings (200 trees, 8 chains, 1000 and 10,000) on three of C2's cells at n 500 (Satellite Image,
  California Housing, Piecewise Linear), 10 replications.
- C4: C1's designs at 200 trees, 25 replications.
- C5: the landed recovery-after-a-response-swap regime, two cells (the large shift at both noise levels), 25
  replications: dbarts' use inside a larger sampler, which no published report measures.

**Readings.** Per cell, as Tan's Figure 14: on C2's four real datasets, RMSE against y and coverage of 95 percent
posterior predictive intervals for y on the test set; on the synthetic designs (C1, C4, C5 and C2's two), RMSE of
f and coverage of 95 percent intervals for f; interval width; R-hat of the test RMSE across chains in Tan's
classical form, with the rank-normalized split R-hat beside it. Beside them, as Pratola, ESS of sigma and of f at
25 test rows, and ESS per second. C5: the landed recovery-sweeps reading. Each as a paired difference with its
interval.

**Reporting.** The report gives every paired difference with its interval and the resolution each cell family
reached, and reads them with no bar, as Tan's L.6 reads its larger contrast. It proposes wording for the news
entry that says what was measured, on which designs, at what resolution; whether to change the entry or the
default is the maintainer's ([For the maintainer](#for-the-maintainer), items 1 and 2).

**Stored results.** The landed study's scripts and per-fit results were not kept; only its tables remain, in
[tree-mixing-proposals.md](../design/tree-mixing-proposals.md). A reread gives nothing beyond those tables: one
chain at 200 trees, 1000 and 2000, 5 seeds, resolving 0.012 on RMSE and 0.02 on coverage.

**Report outline.** 1 page question, the claim and Tan's experiment; 1 page design in ADEMP order; 2 pages readings
(a figure in Tan's form per cell family, a table of paired differences); half a page on C5; half a page on what it
does not settle: one tree and small forests (CGM 1998's regime, where the news entry already says to set swap),
and use inside a larger sampler beyond C5's two cells.

**The study's own choices.** C1 and C5's cells (alternative: Tan's datasets alone); n 500 and 1000 (alternative:
Tan's four sizes, about three times C2's cost with n 10,000); C3 at 10 replications on three cells (alternative:
Tan's full replication, about 10 core-hours); the recovery reading (no published alternative; alternative: drop
C5); the 0.9-34 mixture as the comparator (alternative: add Tan's grow-and-prune arm, about half again C2 and C3's
cost).

## A2 with B4. The logistic k default

**Question.** A user fitting a binary response by logistic BART gets k drawn under chi(1.5, 2). Does that predict
and calibrate probabilities as well as the choices published work makes: CGM's fixed k = 2, Murray's calibrated
fixed scale, and k chosen by cross-validation?

**Published reports.**
- Chipman, George and McCulloch (2010), section 2.2.3 and section 4: k = 2 puts 95 percent prior probability on
  the observed range; on the probit scale, f within (-3, 3).
- CGM (2010), section 5.1 (regression): 42 datasets, 20 random 5/6 splits, BART-default against BART-cv chosen by
  5-fold cross-validation within each training set, relative RMSE per split (an arm's error over the best arm's),
  boxplots and a table of the 50 and 75 percent quantiles.
- CGM (2010), section 5.3, binary: BART-cv over k in {0.25, 0.5, 1, 2, 3} and m in {100, 200, 400, 800}, 20
  replicates, AUC, a one-way ANOVA with replicate as block; "practical differences are not appreciable".
- Murray (2021), section 4.3: a multinomial log-linear prior whose log odds lie within plus or minus 7 with about
  0.95 probability, a0 = 3.5/sqrt(2). Section 5.1, multiclass (3 to 6 categories): 20 UCI datasets, 10 random
  80/20 splits, BART-CV chosen by 10-fold cross-validation within the training set over a0 in {2, 3.5, 6}/sqrt(2),
  accuracy; paired Wilcoxon tests of BART-CV against other methods (Table 1); the default against BART-CV as
  relative accuracy over 2,000 fold comparisons (Figure 4, "within 2 percent ... about 75 percent of the time").
- Hill, Linero and Murray (2020, p. 257): cross-validation improves binary BART, citing Dorie et al. (2019),
  which this plan has not read.

**The translation.** "Murray's k" is the project's own mapping of Murray 4.3 to a normal-leaf binary logit:
Murray's binary log odds are a difference of two log-linear functions, sd sqrt(2) a0 = 3.5, which on dbarts'
logistic scale (k standard deviations from plus or minus pi sqrt(3)) is k = 5.4414/3.5 = 1.555. Murray's grid
translates the same way to k of about 2.72, 1.55 and 0.91.

**Design.**
- Part A, the prior alone, no compute: the induced prior on the logit and on the probability under the default,
  fixed k = 2 and Murray's k, against CGM's and Murray's calibration targets. Run already: under the default the
  95 percent quantile of the absolute logit is 13.4 and 87.5 percent of the prior lies within plus or minus 7;
  fixed k = 2 puts 95 percent within plus or minus 5.33 and 99 percent within 7.
- Part B, the reread, no compute: below.
- Part C, CGM's sections 5.1 and 5.3 with the default as a user gets it: arms BART-default, fixed k = 2, Murray's
  k and BART-cv by `xbart` over CGM's k grid with 5 folds inside each training set, refit at the chosen k. The 22
  datasets the k study read (none on pricing's list) with its 20 splits (10 for the five large ones), 80/20 as
  Murray's (CGM used 5/6), 75 trees. These are the datasets the October study kept the default on, so Part C
  confirms on seen data and is not an out-of-sample test ([For the maintainer](#for-the-maintainer), item 3).
- Part D, calibration on known truths: logistic truths on Friedman and on an additive design, log-odds sd of 0.5,
  1.5 and 3.5, n 200, 1000 and 4000: 18 cells, 100 replications, the same arms less BART-cv.
- Probit, optional: Part D's cells under the probit family, the same arms.

**Readings.** Part C as CGM: raw AUC per arm with CGM 5.3's one-way ANOVA blocked by replicate (Table 5); beside
it, per split, an arm's AUC loss (one minus AUC) over the best arm's, the project's own adaptation of CGM 5.1's
relative form, and the same for log loss and Brier score, which are the study's own choice for "calibrate
probabilities" (alternative: AUC alone, or Murray's accuracy); boxplots over datasets and splits, the table of 50
and 75 percent quantiles; beside them, the default against BART-cv in Murray's Figure 4 form, and a per-dataset
paired Wilcoxon of the default against BART-cv, an adaptation of Murray's Table 1 tests. Log loss near zero on
separable data makes ratios unstable, so the paired differences are shown beside them. Part D: coverage of 95
percent intervals for the probability at test points, RMSE of the probability, bias, width.

**Reporting.** The relative-loss quantiles and Part D's coverage are reported with their intervals and no bar, and
the report recommends as CGM 5.1 reads its default against BART-cv.

**Stored results.** The k study's per-fit summaries hold, for 22 real datasets and 12 simulated cells, log score,
Brier score, coverage, width, probability RMSE and k's R-hat and ESS, for nine priors on k (chi(1, 1),
chi(1.5, 0.5), chi(1.5, 1), the default chi(1.5, 2), chi(1.5, 5), fixed k of 1, 1.5, 2 and 3), at 4 chains of
2000 after 2000, the large datasets subsampled to 4000 rows. No per-row draws. Reread in CGM's form, an arm's log
score over the best of the nine per split, 390 splits: the default's median is 1.0136 and its 75 percent quantile
1.0378; chi(1, 1) is lowest at both, 1.0133 and 1.0369; fixed k = 2 has 1.0300 and 1.2675. On Brier score the
default is lowest at both, 1.0135 and 1.0423. The five drawn priors are not told apart: the default against
chi(1, 1) is better on 8 of 22 datasets, a mean log loss difference of +0.0006 (se 0.0005). All five are far
better than every fixed k in the tail, driven by near-separable data (the default minus fixed k = 2 in log loss:
tictactoe -0.19, sonar -0.082, ionosphere -0.049, banknote -0.037). The stored k chains did not converge (the
default's k R-hat median 1.22, ESS median 13 of 8000 draws) while the predictions largely did (R-hat median
1.005, at most 1.19; ESS median 925). So the reread shows, at no compute, that on these 22 seen datasets at a long
run the drawn priors beat fixed k in the tail; it resolves no ordering among the drawn priors, and a best-of-nine
ratio depends on which arms are in the set. It does not settle the default run length, BART-cv, Murray's k or
calibration on known truths.

**Report outline.** 1 page question and published calibrations (Part A's figure); 1 page Part B and C's
boxplots and quantile table; 1.5 pages Part D's coverage figure; 1 page what it means and does not settle: that
the drawn default beats other drawn priors, k's mixing, and probit.

**The study's own choices.** The translation of Murray's prior; the best-of-nine reread (alternative: CGM's
best-of-arms over Part C's arms only); m not cross-validated, held at 75 (alternative: CGM's m grid, about 4
times Part C); Part C run before pricing's ruling, since its readings are predictions, which converge sooner
than k (alternative: wait); Part D's truths and sizes (alternative: CGM's section 6 simulation, which is gaussian
only); probit left optional (alternative: run it, about 4 more core-hours at Part D's size).

## A1. The negative binomial default

**Question.** A user fitting counts with `family = nbinom` gets the default leaf prior on the log mean. On count
data people meet, does it estimate the means, group means, rate ratios and new counts they report, with RMSE,
bias and coverage as good as Murray's calibration?

**Published reports.**
- Murray (2021), section 4.4: the leaf prior on the log mean calibrated to the data, a0 = 0.5 [log(y*) - log(mu0)]
  for "a near-maximum value (or upper quantile of the empirical distribution) of y" y* and a base mean mu0, and
  "for large values of mu0 it may be necessary to inflate this value"; section 5.2, an example with a0 = 2 by
  "slightly inflating" the heuristic at y* = 50; models compared by WAIC and 10-fold cross-validated log loss
  (Table 3), with no cross-validated choice of a0.
- Hill (2011): real covariates, outcomes simulated from response surfaces written down (random coefficients), 1000
  replications.
- The referees of the landed study asked for natural estimands per dataset, bias, and more datasets chosen by a
  stated rule.

**Murray's calibration, as read here** (the study's own reading of 4.4, fixed before any fit): mu0 is the
default's own center, the total count over the total exposure; y* is the largest count per unit exposure in the
training rows; a0 = 0.5 log(y*/mu0), not inflated; the arm's k is 3/a0, dbarts' k.scale of 3 on the log mean.
Alternatives: y* at the 97.5 percent quantile; Murray's inflation to the next half unit, as in 5.2.

**Design.**
- Part A, the reread: below.
- Part B, no compute: Murray's a0 per real dataset, and where it falls under the default's induced prior on the
  spread of the log mean.
- Estimands, before any fit: for each dataset, a table of its source, the publication its documentation cites,
  and the estimand that publication reports (a rate ratio for a named treatment or exposure, means of named
  groups, or the fitted mean where it reports neither), read from the publication itself and checked by the
  harness review before any fit. The table is not written here: this plan's author has not read those
  publications.
- Part C: arms the default, Murray's calibration fixed per dataset and a drawn prior centered at Murray's value;
  the 11 real datasets with the landed study's two truths each, 20 replications. Real outcomes by Murray's
  10-fold log loss and WAIC, the same three arms ([For the maintainer](#for-the-maintainer), item 5, on a
  cross-validated arm).
- The two truths are not neutral, and are kept as the landed study's: the BART truth is fitted at fixed k = 1,
  rougher than the default's prior center (median k 1.9); the as-fitted truth is a `glm.nb` or `gam` additive fit,
  smoother. Both are fitted to each dataset's outcome, the study's own form; the published form surveyed is
  Hill's, surfaces written down on real covariates. Coverage differences between arms can then partly reflect
  the generator, so every table is per dataset and per generator, never pooled.
- If the maintainer adds datasets, the added ones run Part C at the same size.

**Readings.** Per dataset, generator and estimand: RMSE, bias, coverage of 90 percent intervals and width, with
intervals, on its own spread; real outcomes as Murray's Table 3, log loss and WAIC per dataset.

**Reporting.** The default and Murray's calibration side by side per dataset and generator, with paired intervals
and no bar; the report names the datasets where either is ahead beyond its interval, and recommends.

**Stored results.** The landed study's scores hold, for the default at its default length, 5 replications per
real dataset (10 for crimtab), coverage of 50, 80, 90 and 95 percent intervals, interval score, width, RMSE
and log density by estimand, for each truth and its widened, subsampled and fold cells. Its full draws, with
truths, are archived on the laptop, so bias reads at no compute. Reread of fitted means: the default's 90 percent
coverage runs from 0.762 (Traffic, BART-fitted truth) to 1.000 (Insurance), and of rate ratios from 0.789 to
1.000. Between-replication standard errors are 0 to 0.033, so the reread resolves cells: four cells cover below
0.90 beyond a t(4) interval, all with the BART-fitted truth (Traffic 0.762 [0.669, 0.855], Arabidopsis 0.820
[0.750, 0.890], mite 0.852 [0.834, 0.871], solder 0.853 [0.822, 0.885]), and 13 cover above 0.90 beyond it. Every
as-fitted truth covers above 0.90 at its mean (roaches lowest, 0.911), and 7 of 11 BART-fitted truths below it.
The long runs' same cells move some and not others: Traffic 0.778 and solder 0.866 stay below 0.90 beyond their
intervals, mite reaches 0.895. The reread gives every dataset's own table with bias added; it reads the default
alone, so it cannot say whether Murray's calibration does better.

**Report outline.** 1 page question, Murray's calibration and Part B's figure; 2 pages per-dataset tables
(RMSE, bias, coverage, width) by estimand and generator; 1 page real outcomes in Murray's form; 1 page meaning and
limits: no claim over count data at large from 11 datasets, and none that attributes coverage to the prior where
the generators differ.

**The study's own choices.** The reading of Murray's y* and mu0; the drawn prior centered at Murray (alternative:
Murray's fixed arm alone); 20 replications (alternative: 50, about 2.5 times Part C's cost); the two fitted
truths (alternatives: as fitted only; Hill-style surfaces written on each dataset's covariates).

## A3 with A4. stan4bart's initial fit and the starting sigma on every forest's basis

### A3. The initial fit

**Question.** A stan4bart user gets the forest's range from an initial linear mixed model that holds the forest's
columns as linear terms. Does that, against splines on those columns, change the treatment effects, group effects
and intervals they report, across designs and on real multilevel data?

**Published reports.**
- Dorie, Perrett, Hill and Goodrich (2022), sections 6 and 7: IHDP covariates, grouping by site (8) or by
  mother's age (26), varying intercepts or intercepts and slopes, response surfaces A, B and C as Hill (2011): 12
  settings. Effect sizes 0, 0.2, 0.5 and 0.8 outcome sds and intraclass correlations 0.2, 0.333 and 0.5 were
  explored and not shown: "We found no difference in results across these choices" (6.2.2). Readings: SATT,
  group ATT and PEHE, RMSE and interval length divided by the outcome's sd, 95 percent coverage. The replication
  count is not stated in the text read.
- Hill (2011): 1000 replications of IHDP surfaces.

**Design.** Arms: linear (shipped) and three-degree-of-freedom splines, the landed study's two. Cells: the paper's
12 IHDP settings at ICC 0.333 and effect size 0.5, one of each on the paper's own finding of no difference across
them, 200 replications; the landed study's half-curved family, 100 seeds; radon (919 homes, 85 counties) and High
School and Beyond (7185 students, 160 schools), 200 and 50 replications, each with two forest-fitted truths, one
from a stan4bart fit at each arm's initial fit (default length, posterior mean), so that neither arm generates
every truth; both datasets' real outcomes by 10-fold cross-validation.

**Readings.** The paper's: SATT and group ATT RMSE over the outcome sd, PEHE, coverage and interval length; on
radon and HSB, per generator, county or school means and the fixed effects (log uranium; sector, mean SES), RMSE,
bias and coverage; real outcomes by held-out log density.

**Reporting.** Paired differences between the arms per cell with their intervals and no bar, read as the paper
reads its methods across its 12 settings; the report recommends linear, splines or a documented knob.

**Stored results.** The landed study's 240 in-between fits and its earlier stress fits are not in their summary
form; its tables remain. No reread beyond them.

### A4. The starting sigma on every forest's basis

**Question.** A user writing a two-forest treatment model gets a starting sigma from a regression on x and the
treatment basis. Does that, against the regression on x alone and against the sample sd, change the effect
estimates, their intervals and the noise sd's coverage?

**Published reports.**
- CGM (2010), section 2.2.4: the sigma prior is set from "a rough data-based overestimate": the naive sample sd
  of y, or the residual sd of a linear regression; q = 0.90 by default. The estimate calibrates the prior, as here.
- The bcf package: the residual sd of `lm(y ~ z + x)`, nu 3, q 0.9.
- Hahn, Murray and Carvalho (2020), section 6.1: 8 designs (n 250 or 500, homogeneous or heterogeneous effect,
  linear or nonlinear prognostic function), 200 replications, ATE and CATE RMSE, coverage and interval length
  (Tables 2 and 3). The section does not state the noise sd.

**Design.** Arms: shipped (x and the basis, bcf's choice), x alone (the old default, by an explicit `sigest`) and
the sample sd (CGM's naive choice). Cells: BCF's 8 designs with noise sd 1, the study's own value (alternative:
noise sd 0.5 and 2 as a factor, three times the cost), 200 replications; IHDP surface B at 200. Where 6.1's text
names the three-level variable x4 in g(x4) while defining x5 as the categorical one, g is applied to the
three-level variable (the study's own reading).

**Readings.** BCF's Tables 2 and 3: ATE and CATE RMSE, coverage and length of their intervals. Beside them, two
readings of the study's own, which BCF reads neither of: coverage of the noise sd's interval, and the share of
replications where the starting estimate is at least the true sd (CGM's intent of an overestimate).

**Reporting.** Each reading per arm and the paired differences, with intervals and no bar and no expected result;
the report recommends.

**Stored results.** The landed 20-seed check is summarized in its design note only. No reread.

**Report outline (A3 and A4 together).** 1 page questions and published designs; 1.5 pages A3 in the paper's
tables; 1 page A4 in BCF's tables; 1 page radon and HSB; half a page limits: the paper's 12 settings at one ICC
and effect size; nothing about stan4bart 0.0-13's initial fit unless that arm loads; nothing about the starting
sigma in other families or bases.

**The study's own choices.** A3: one ICC and effect size (alternative: all three ICCs, about 3 times the IHDP
cost); 200 replications (alternative: Hill's 1000); forest-fitted truths for radon and HSB, one per arm
(alternatives: real outcomes only; a truth from a third model). A4: the noise sd of 1; the two readings BCF does
not read.

## A6. The monotone defaults

**Question.** A user constraining f to be monotone gets the "joint" prior at cgm()'s tree prior. Do 95 percent
intervals for f cover near nominal, and is the fit's tree size near what an unconstrained fit gives, which mBART
aims at for its own approximation?

**Published reports.**
- Chipman, George, McCulloch and Shively (2022), mBART. Its stated model keeps CGM's tree prior, p(T*)/p(T0)
  "the same as in the unconstrained case" (4.2), with a normalized constrained leaf density (4.13) and the
  constrained leaves' variance inflated by c^2 = pi/(pi - 1) (3.3). Its sampler omits the normalizer and
  compensates with alpha .25 and beta .8, reaching "tree sizes comparable to those obtained in unconstrained BART"
  (4.3).
- mBART section 5.3: f = x1 x2^2 + x3 x4^3 + x5 on five uniform predictors (read from the rendered arXiv page;
  the Bayesian Analysis version, DOI 10.1214/21-BA1259, is checked before coding), sigma in {0.2, 0.5, 0.7, 1, 2},
  200 datasets of 500 training and 1000 test rows; RMSE ratio boxplots (Figure 8) and Table 1's width and
  coverage; mBART covered 90.1 to 97.9 percent, read as "good reliability", with no bar.

**Under an exact sampler.** dbarts applies mBART's c^2 to every constrained leaf under either prior (read in the
engine's source). Its "leaf" prior at cgm() defaults is then the nearest counterpart of mBART's stated model, and
"joint" at (.25, .8) of what mBART's software samples. The shipped default is "joint" at cgm() defaults
(dec-B151, dec-B152).

**Design.** Arms: unconstrained; "joint" at cgm() defaults (shipped); "joint" at mBART's tree prior (.25, .8),
which is mBART's tree prior only. Whether "leaf" joins is the maintainer's
([For the maintainer](#for-the-maintainer), item 7).
- C1, the reread: below.
- C2, the partial rerun, at 75 trees and 4 chains: the three low-noise cells at n 2000, rerun because every arm,
  the unconstrained included, covers below nominal there at 50 trees on one chain of 500, so the default's 4
  chains are the question; 50 replications; and mBART's section 5.3 at its five noise levels, 50 replications.

**Readings.** As mBART's Table 1 and Figure 8: coverage and width of 95 percent intervals for f, RMSE ratios
against the unconstrained fit, and splits per tree against the unconstrained fit; with intervals.

**Reporting.** Each arm's readings and the paired differences, with intervals and no bar; the unconstrained fit
covers 0.96 to 0.99 in 16 of the 24 stored cells, so coverage is read against nominal, not against the other arms.

**Stored results.** The tree-prior study's baseline holds 24 cells (three designs, two noise levels, n 200 and
2000, 50 and 200 trees), 8 replications, one chain of 500 after 500. Reread: in all 24 cells coverage orders
unconstrained, then shipped, then mBART's tree prior; by paired t(7) intervals the shipped prior covers above
mBART's in 23 of 24 cells and the unconstrained fit above the shipped one in 16. Splits per tree, per cell: shipped
1.14 to 1.47, unconstrained 1.23 to 1.56, mBART's tree prior 0.25 to 0.96 (averaged over the six cells of each
noise and tree count: 1.27 to 1.32, 1.40 to 1.45 and 0.29 to 0.68). The shipped prior is nearer the unconstrained
fit's tree size than mBART's tree prior is in every cell. At 50 trees, low noise and n 2000 every arm covers below
nominal (the unconstrained fit 0.835 to 0.869), which is not specific to the constraint; one chain of 500 is also
a candidate cause. Not settled: the default's 4 chains and 75 trees, mBART's design, and "joint" against "leaf".

**Report outline.** 1 page question and mBART's calibration; 1.5 pages the reread's table and tree-size figure;
1.5 pages C2 in mBART's Table 1 and Figure 8 form; 1 page meaning and limits.

**The study's own choices.** 50 replications on mBART's design (alternative: its 200, about 4 times C2's cost);
the rerun cells (alternative: all 12 stored cells at 75 trees).

## Sizing

Costs are core-hours on the box. Half-widths are of a 95 percent interval at first look and at full size: by eq.
(1) for single-interval coverage; by Table 6 with s, the per-replication spread from stored runs, for coverage
over test points ([Common to every study](#common-to-every-study)). Targets are the project's own.

| study | component | size (looks) | expected | at most | half-width, first / full | basis | target | full run needs |
|---|---|---|---|---|---|---|---|---|
| A5 | C1 | 5 cells x 2 arms x 25 (13, 25) | 0.2 | 0.3 | 0.009 / 0.006 | s 0.016, landed 5 seeds | 0.02 | met at 13 |
| A5 | C2 | 12 cells x 2 x 25 (13, 25) | 0.3 | 0.5 | from smoke | none stored | 0.02 | met at 25 if s is at most 0.05 |
| A5 | C3 | 3 cells x 2 x 10 (5, 10) | 2.0 | 3.0 | R-hat, reported | | none | full size |
| A5 | C4 | 5 x 2 x 25 (13, 25) | 0.6 | 0.9 | as C1 | as C1 | 0.02 | as C1 |
| A5 | C5 | 2 x 2 x 25 (13, 25) | 0.2 | 0.3 | 10% of sweeps at full | landed spread | 20% | met at 25 |
| A2 | A, B | none | 0 | 0 | | | | |
| A2 | C | 22 datasets x 20 splits (10, 20) | 3.5 | 5.3 | median relative loss 0.005 at full | stored splits | 0.01 | met at 20 |
| A2 | D | 18 cells x 3 x 100 (50, 100) | 7.0 | 10.5 | 0.028 / 0.020 at median s 0.10; 0.053 / 0.037 at s 0.19 | k study's simulated cells | 0.02 | 350 at s 0.19 |
| A1 | A, B | none | 0.1 | 0.1 | | | | |
| A1 | C, truths | 22 cells x 3 x 20 (10, 20) | 2.0 | 3.0 | fitted means 0.016 / 0.011 at s 0.026; group means 0.068 / 0.048 at s 0.11 | stored replications | 0.05 | 148 at the largest s, the rate ratios' 0.31 |
| A1 | C, folds | 11 x 10 folds x 3 arms | 1.3 | 1.9 | log loss, reported | | none | full size |
| A3 | IHDP | 12 x 2 x 200 (100, 200) | 4.0 | 6.0 | SATT 0.043 / 0.030 | eq. (1) | 0.03 | met at 200 |
| A3 | half-curved | 100 seeds x 2 (50, 100) | 0.3 | 0.5 | group-effect error 2% | landed spread | 3% | met at 100 |
| A3 | radon | 200 x 2 arms x 2 truths (100, 200) | 0.6 | 0.9 | fixed effect 0.043 / 0.030 | eq. (1) | 0.03 | met at 200 |
| A3 | HSB | 50 x 2 x 2 (25, 50) | 1.1 | 1.7 | fixed effect 0.085 / 0.060 | eq. (1) | 0.03 | 200, about 4.4 core-hours |
| A3 | folds | radon and HSB, 10 folds x 2 | 0.2 | 0.3 | log density, reported | | none | full size |
| A4 | BCF | 8 x 3 x 200 (100, 200) | 1.9 | 2.9 | ATE and noise sd 0.043 / 0.030 | eq. (1) | 0.03 | met at 200 |
| A4 | IHDP B | 3 x 200 (100, 200) | 0.6 | 0.9 | as BCF | eq. (1) | 0.03 | met at 200 |
| A6 | C1 | reread | 0 | 0 | | | | |
| A6 | C2 | 8 cells x 3 x 50 (25, 50) | 5.0 | 7.5 | 0.005 / 0.004 at median s 0.013; 0.024 / 0.017 at s 0.06 | stored baseline | 0.02 | met at 50 |

A first look cannot meet its target in A2's Part D, A1's group means, A3's IHDP and radon, A4 and A6's widest
cells; the last column says what full size needs. A4's paired difference in ATE coverage resolves 0.031 at 200
replications if the arms disagree on at most 5 percent of them, 0.044 at 10 percent.

Expected 30.9 core-hours (A5 3.3, A2 10.5, A1 3.4, A3 6.2, A4 2.5, A6 5.0), at most 46.5. The options for the
maintainer add: "leaf" in A6, about 1.7 (2.5 at most); a nested cross-validated arm in A1, about 8 (12); A2's
fresh datasets about 2 and its wider CV grid about 1.5; HSB at 200 about 3.3; added datasets for A1 about 0.5
each. Timings used: a default gaussian fit 1.8 to 2.7 s at n 500 to 2000, probit 1.6 s at n 2000, logistic about
2, 4 and 8 s at n 500, 2000 and 4000, negative binomial a median of 4 s (2 to 25 s by dataset), stan4bart on radon
2.6 s (HSB taken at 8 times radon, its rows' ratio, to be timed in the smoke), a two-forest treatment fit 1.2 s at
n 500, a constrained monotone sweep 4 to 5 ms at 50 trees, "leaf" 1.04 to 1.06 times "joint" per sweep at 50 and
200 trees (the feel study).

## For the maintainer

Each changes what a study measures, its comparators, or a choice the maintainer ruled; the plan is written to the
recommendation and nothing runs on it until ruled.

1. A5, C2's readings on real data. (a) Tan's L.1: RMSE against y and coverage of posterior predictive intervals
   for y on the four real datasets, f on the synthetic ones; recommended, as written; no cost. (b) Outcomes
   simulated from surfaces written on Tan's real covariates, Hill's form, reading f everywhere; about 0.1 more
   core-hours and a generator of the study's own.
2. A5, the news entry. "Costs nothing measurable" names no resolution. (a) Reword it after the run to what was
   measured, on which designs, at what resolution, whatever the readings; recommended. (b) Keep it unless a
   reading moves. A change of default, if the report recommends one, is also yours.
3. A2, Part C's datasets. (a) The 22 seen datasets alone, reported as a confirmation on seen data; no cost. (b)
   Add about 12 fresh binary datasets drawn by the rule independent agents are building for pricing under
   dec-B635, once ruled, none from pricing's held-out list; about 2 core-hours; recommended, since a default kept on
   these 22 datasets is then tested on others.
4. A2, the BART-cv grid. (a) CGM's k grid {0.25, 0.5, 1, 2, 3}, as written. (b) Murray's grid translated,
   {0.91, 1.55, 2.72}. (c) Both, 8 values; about 1.5 more core-hours; recommended, as each is a published grid.
5. A1, a cross-validated arm. (a) None, as Murray scores counts, with the default against Murray's calibration;
   recommended, as written. (b) k chosen over {0.5, 1, 2, 4} by 5-fold cross-validation inside each of the 10
   outer training folds, as CGM 5.1 and Murray 5.1 select, refit at the choice; about 8 more core-hours, 12 at most.
   Choosing k on the outer folds themselves would report the best k in hindsight, and is not offered.
6. A1's added datasets. (a) Independent agents choose about 12 count datasets by a stated rule, as dec-B635 does
   for pricing; (b) the maintainer names them; (c) the 11 alone. Recommend (a): the referees asked for it, and a
   rule is what published comparisons state. About 0.5 core-hours each.
7. A6, the "leaf" prior. The study is titled the monotone defaults, and "joint" is one of them (dec-B151, ruled
   on the feel study's one small design). (a) Add "leaf" at cgm() defaults to C2, the nearest counterpart of
   mBART's stated model; about 1.7 core-hours; recommended; its readings bear on dec-B151. (b) Narrow the study's
   title and claims to the tree prior under "joint"; no cost.
8. A3, stan4bart 0.0-13 as an arm. It does not load beside today's libraries; a private build against dbarts
   0.9-34 on the box might. Recommend one attempt, the arm added if it loads and omitted, with that said, if not.
9. A3, HSB's replications. At 50, HSB's fixed-effect coverage reaches a half-width of 0.060 against its 0.03
   target. (a) 50, as written, the width reached reported; no cost. (b) 200, which meets the target; about 3.3
   more core-hours (HSB's fit time a guess until the smoke). Recommend (b) if the budget in item 10 allows it.
10. Run order and budget. dec-B637's ruled first tier was shown as "the cheap items now, about 5 core-hours",
    larger reruns brought one at a time; this plan asks about 31 core-hours expected and 46.5 at most as written,
    roughly six times that, and about 39 and 59 with items 3, 4, 7 and 9, before A1's added datasets. Recommend
    the order above, A2's Part D after pricing's ruling, and approval one study at a time against that figure.

Note: where a report recommends a change of default or of a documented claim (A5's news entry, A2's k, A1's
calibration, A6's prior), that change is the maintainer's; each report recommends one.

## Referee round 1

One referee accepted all six studies with changes. Each required change, and where it is answered:

- Common: 1, every "What results mean" is now a Reporting paragraph with no rule; 2, Morris et al. for sizing and
  Monte Carlo SEs (Common; [Sizing](#sizing)); 3, single-interval coverage resized by eq. (1), the paired
  difference with its discordance. Optional 4 (no adjusted intervals) and 5 (ADEMP order) taken.
- A5: 1, Tan's contrast named the nearest and a larger one; 2, y on real data, f on synthetic (item 1); 3, the
  rule dropped, the resolution reached reported, the news wording item 2. Optional 4 and 5 taken.
- A2: 1, the reread reworded, with k's and the predictions' convergence; 2, the translation named the project's
  own, Murray's grid named (item 4); 3, the Wilcoxon named an adaptation, Figure 4 and multiclass stated; 4, Part
  C named a confirmation on seen data (item 3); 5, log loss and Brier named the study's own, AUC beside them.
  Optional 6 and 7 taken.
- A1: 1, the generators stated and tabled apart; Hill (2011) writes its surfaces on real covariates rather than
  fitting them to outcomes, so the fitted truths are named the study's own; 2, the cross-validated arm dropped,
  nested selection offered (item 5); 3, y* and mu0 fixed; 4, the estimand table a step before any fit, and why it
  is not written here; 5, rule removed. Optional 6 taken (the long runs' cells).
- A3: 1, IHDP at 200; 2, effect size 0.5, on the paper's own finding; 3, one truth from each arm; 4, rule removed.
- A4: 1, noise sd 1 named; 2, 200 replications for every arm; 3, no expected result; 4, both readings named.
- A6: 1, "leaf" is item 7, the claim narrowed if declined; 2, the mBART arm named its tree prior only. The
  referee found no c^2 in the package, but the engine applies it to every constrained leaf under either prior;
  3, splits ranges and formula corrected; 4, the settle rule dropped, the rerun cells' reason stated; 5, both
  claims withdrawn; 6, rule removed.
- Numbers, by rerunning the referee's reread scripts: the k reread (chi(1, 1) has the lowest median, 1.0133,
  against the default's 1.0136), the negative binomial reread (the cells it resolves), and the monotone reread
  (per-cell splits ranges; 23 and 16 of 24 cells by paired intervals).

## Papers

Read in full: Chipman, George and McCulloch (1998, 2010); Chipman, George, McCulloch and Shively (2022); Tan,
Ronen, Saarinen and Yu (2024); Murray (2021, the arXiv version); Pratola (2016); He and Hahn (2023); Hahn, Murray
and Carvalho (2020); Hill (2011); Hill, Linero and Murray (2020); Vehtari et al. (2021); Dorie, Perrett, Hill and
Goodrich (2022, the full text from Europe PMC; the publisher's PDF was refused); Morris, White and Crowther (2019,
arXiv 1712.03198; DOI 10.1002/sim.8086). Not reached: mBART's used-car data (section 5.4), whose availability is
unverified; the stan4bart paper's replication count, not in the text; Dorie et al. (2019), cited through Hill,
Linero and Murray; the Bayesian Analysis version of mBART.

## Verification

- Run: 0.9-34's tree-move mixture and binary k prior on the released package; the induced logit prior under
  the default, fixed k = 2 and Murray's k; the rereads of the k study's real splits, of the negative binomial
  scores (default length and long runs) and of the monotone baseline quoted above; the per-replication spreads in
  [Sizing](#sizing); the timings there.
- Read: the bcf package's starting sigma; stan4bart 0.0-13's initial fit, in its source; the engine's constrained
  leaf variance.

## Not verified

- Box timings: every cost is from the laptop under load, to be rechecked in each smoke; HSB's fit time is a guess.
- That PMLB's copies of Tan's four datasets match Tan's, and that Tan's synthetic designs are copied exactly.
- The spreads behind the sizing: from stored runs at other settings (A2's from probit truths, A5's from 5 seeds
  at 200 trees), rechecked at each smoke.
- That the negative binomial scores reread are the rescoring after the rate-ratio fix; that fix touched designed
  cells, and only real-data cells are quoted.
- That the archived negative binomial draws read back whole.
