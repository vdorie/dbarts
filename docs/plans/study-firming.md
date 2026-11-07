# study-firming: the first tier of landed studies brought to published practice

Status: planned (2026-11-07). Nothing has run. Five studies, each separable, in the order the maintainer set
(dec-B637): the swap move, the logistic k default, the negative binomial default, stan4bart's initial fit with the
starting sigma, and the monotone defaults. Each study's design, readings and rules are written here before any fit
it reads. A rule says what a reading shows and recommends; the maintainer decides from a report of at most 5 pages
of text (dec-B559).

agent: per study, one opus implementer for the harness and the runs on the bench box, one review of the harness
and its smoke (the gate-only tier of [Process by risk](README.md#process-by-risk)), and two independent referees
on the report. The rereads need no implementer: one agent on the laptop.
rng: the studies change no package code. Every run is on one 1.0 library in shipped mode; an arm is a user-level
setting, never a build.
window: before 1.0-0. The logistic k study's Part D waits for k-move-pricing's ruling on the move (dec-B634).
Studies run one at a time on the box.
budget: about 23 core-hours expected over the five, about 35 at most, 29 and 44 with A1's added datasets
([Sizing](#sizing)); each study's harness a few hundred lines of R, much of it reused from the landed studies'.

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
- Paired arms: every arm of a replication shares its data seed and its sampler seed.
- Every reading is a mean over replications with a 95 percent interval; a difference between arms is a paired
  difference with its 95 percent interval. Where a study reads many cells, Holm-adjusted intervals are reported
  beside the plain ones (the project's own; the alternative is plain intervals alone, as every report surveyed
  shows them).
- Machine: the x86 bench box for every run, the laptop for every reread. Costs below are from single-thread
  timings on the laptop under a load of 8 to 11 and are rechecked in each harness's smoke.
- Stop rule (the project's own, as pricing's is): two looks, at half and full size. A component stops at the
  first look once every interval it reports is narrower than its resolution target in the sizing table; the
  targets are the project's own, set to the size of difference the landed study could resolve. No component runs
  past 1.5 times its expected cost: it stops and reports what it has.
- Each report: at most 5 pages of text; each table or figure tied to the published table or figure it copies;
  the study's own choices listed with their alternatives; a section on what it does not settle.

## A5. The swap move dropped from the default mixture

**Question.** A user who never sets `proposal.probs` gets a mixture without the swap move. Does that change the
accuracy, the interval coverage or the mixing of what they estimate, against the 0.9-34 mixture? The news entry
says "dropping it costs nothing measurable in accuracy, coverage or mixing".

**Published reports.**
- Tan, Ronen, Saarinen and Yu (2024), appendix L.6, "Experiment 7: Restricted vs. full move set": the same
  contrast, the full mixture (grow .25, prune .25, change .4, swap .1) against grow and prune alone. 200 trees,
  8 chains, 1000 burn-in and 10,000 kept, 25 replicates; training sizes 200, 500, 1000 and 10,000 drawn from
  four PMLB datasets (Breast Tumor, Echo Months, Satellite Image, California Housing; a fixed 10 percent test
  set) and two synthetic designs (Low-Dimensional Smooth, Piecewise Linear). Readings: Gelman-Rubin R-hat of
  the test RMSE across chains, 95 percent coverage and RMSE, with bars of 1.96 standard errors. Verdict, read
  with no bar: "restricting the move set does not substantially affect R-hat, coverage, or RMSE."
- Chipman, George and McCulloch (1998), single tree: "omission of the SWAP step inhibits mixing".
- Pratola (2016): effective sample size of sigma and per second (Table 2), run times (Table 3).
- He and Hahn (2023), Table 4: coverage, interval length, RMSE and time over 100 replications.
- Vehtari et al. (2021): 4 chains, R-hat below 1.01 and bulk ESS of at least 400 as the convergence check.

**Design.** Two arms: the shipped mixture and the 0.9-34 mixture, both set through `dbartsControl` on one library.
- C1: the landed study's five designs (Friedman at two noise levels, a depth-3 interaction, 100 predictors, a
  probit), 25 replications at the defaults.
- C2: Tan's six datasets at n 500 and 2000, 25 replications at the defaults, Tan's test sets.
- C3: Tan's settings (200 trees, 8 chains, 1000 and 10,000) on three of C2's cells at n 500 (Satellite Image,
  California Housing, Piecewise Linear), 10 replications.
- C4: C1's designs at 200 trees, 25 replications.
- C5: the landed recovery-after-a-response-swap regime, two cells (the large shift at both noise levels), 25
  replications: dbarts' use inside a larger sampler, which no published report measures.

**Readings.** Per cell, as Tan's Figure 14: R-hat of the test RMSE across chains, coverage of 95 percent
intervals for f on the test set, RMSE of f, interval width; beside them, as Pratola, ESS of sigma and of f at
25 test rows, and ESS per second. C5: the landed recovery-sweeps reading. Each as a paired difference with its
interval.

**What results mean.** Every paired interval containing zero, or every interval favoring the shipped mixture:
the news claim stands as written. An interval excluding zero against the shipped mixture in some cell: the
report names the cell, and the recommendation is to reword the claim to what was measured; a change of default
is recommended only if the same cell's coverage or RMSE moves by more than the landed study could resolve.

**Stored results.** The landed study's scripts and per-fit results were not kept; only its tables remain, in
[tree-mixing-proposals.md](../design/tree-mixing-proposals.md). A reread gives nothing beyond those tables: one
chain at 200 trees, 1000 and 2000, 5 seeds, resolving 0.012 on RMSE and 0.02 on coverage.

**Report outline.** 1 page question, the claim and Tan's experiment; 1 page design; 2 pages readings (a figure
in Tan's form per cell family, a table of paired differences); half a page on C5; half a page on what it does
not settle.

**The study's own choices.** C1 and C5's cells (alternative: Tan's datasets alone); n 500 and 2000 (alternative:
Tan's four sizes, about three times C2's cost with n 10,000); C3 at 10 replications on three cells (alternative:
Tan's full replication, about 10 core-hours); the recovery reading (no published alternative; alternative: drop
C5).

## A2 with B4. The logistic k default

**Question.** A user fitting a binary response by logistic BART gets k drawn under chi(1.5, 2). Does that predict
and calibrate probabilities as well as the choices published work makes: CGM's fixed k = 2, Murray's calibrated
fixed scale, and k chosen by cross-validation?

**Published reports.**
- Chipman, George and McCulloch (2010), section 2.2.3 and section 4: k = 2 puts 95 percent prior probability on
  the observed range; on the probit scale, f within (-3, 3).
- CGM (2010), section 5.1: 42 datasets, 20 random 5/6 splits, BART-default against BART-cv chosen by 5-fold
  cross-validation, relative RMSE per split (an arm's error over the best arm's), boxplots and a table of the 50
  and 75 percent quantiles.
- CGM (2010), section 5.3, binary: BART-cv over k in {0.25, 0.5, 1, 2, 3} and m in {100, 200, 400, 800}, 20
  replicates, AUC, a one-way ANOVA with replicate as block; "practical differences are not appreciable".
- Murray (2021), section 4.3: the leaf prior calibrated so log odds lie within plus or minus 7, a0 = 3.5/sqrt(2),
  which is k of about 1.555 on the logistic scale; section 5.1: 20 UCI datasets, 10 random 80/20 splits, 10-fold
  cross-validation, accuracy, paired Wilcoxon tests and relative-accuracy figures ("within 2 percent ... about 75
  percent of the time").
- Hill, Linero and Murray (2020): cross-validation improves binary BART.

**Design.**
- Part A, the prior alone, no compute: the induced prior on the logit and on the probability under the default,
  fixed k = 2 and Murray's k, against CGM's and Murray's calibration targets. Run already: under the default the
  95 percent quantile of the absolute logit is 13.4 and 87 percent of the prior lies within plus or minus 7;
  fixed k = 2 puts 95 percent within plus or minus 5.3 and 99 percent within 7.
- Part B, the reread, no compute: below.
- Part C, CGM's section 5.1 and 5.3 with the default as a user gets it: arms BART-default, fixed k = 2, Murray's
  k and BART-cv by `xbart` over k in {0.25, 0.5, 1, 2, 3} with 5 folds, refit at the chosen k. The 22 datasets
  the k study already read (none on pricing's list), its 20 splits (10 for the five large ones), 75 trees.
- Part D, calibration on known truths: logistic truths on Friedman and on an additive design, log-odds sd of 0.5,
  1.5 and 3.5, n 200, 1000 and 4000: 18 cells, 50 replications, the same arms less BART-cv.
- Probit, optional: Part D's cells under the probit family, the same arms.

**Readings.** Part C as CGM: per split, an arm's log loss over the best arm's on that split, and the same for
Brier score and AUC; boxplots over datasets and splits, the table of 50 and 75 percent quantiles; beside them,
Murray's per-dataset paired Wilcoxon of the default against BART-cv, and CGM's ANOVA. Log loss near zero on
separable data makes ratios unstable, so the paired differences are shown beside them. Part D: coverage of 95
percent intervals for the probability at test points, RMSE of the probability, bias, width; no bar.

**What results mean.** The default's relative-loss quantiles within those of BART-cv and fixed k = 2 on Part C,
and its coverage in Part D no further from 0.95 than the fixed arms': the default stands as CGM's comparison would
keep it. A fixed arm better on both: the recommendation names it, and the maintainer decides.

**Stored results.** The k study's per-fit summaries hold, for 22 real datasets and 12 simulated cells, log score,
Brier score, coverage, width, probability RMSE and k's R-hat and ESS, for nine priors on k (chi(1, 1),
chi(1.5, 0.5), chi(1.5, 1), the default chi(1.5, 2), chi(1.5, 5), fixed k of 1, 1.5, 2 and 3), at 4 chains of
2000 after 2000, the large datasets subsampled to 4000 rows. No per-row draws. Reread in CGM's form, the arm's
log score over the best of the nine per split: the default's median is 1.014 and its 75 percent quantile 1.038,
the lowest median of the nine and within 0.001 of the lowest at 75 percent; fixed k = 2 has 1.030 and 1.268, with
a long tail from near-separable data. Brier score agrees. That settles, at no compute, that the drawn default is
at least as good as any of the nine on real outcomes at a long run; it does not settle the default run length,
BART-cv, Murray's k or calibration on known truths.

**Report outline.** 1 page question and published calibrations (Part A's figure); 1 page Part B and C's
boxplots and quantile table; 1.5 pages Part D's coverage figure; 1 page what it means and does not settle.

**The study's own choices.** The best-of-nine reread (alternative: CGM's best-of-arms over Part C's arms
only); m not cross-validated, held at 75 (alternative: CGM's m grid, about 4 times Part C); Part C run before
pricing's ruling, since its readings are predictions, which converge sooner than k (alternative: wait); Part D's
truths and sizes (alternative: CGM's section 6 simulation, which is gaussian only); probit left optional
(alternative: run it, about 2 more core-hours).

## A1. The negative binomial default

**Question.** A user fitting counts with `family = nbinom` gets the default leaf prior on the log mean. On count
data people meet, does it estimate the means, group means, rate ratios and new counts they report, with RMSE,
bias and coverage as good as Murray's calibration and a cross-validated choice?

**Published reports.**
- Murray (2021), section 4.4: the negative binomial leaf prior calibrated to the data, a0 = 0.5 (log y* - log
  mu0) for a high count y* and a base mean mu0; section 5.2 (an example with a0 = 2): models compared by WAIC and
  10-fold cross-validated log loss (Table 3).
- The referees of the landed study asked for natural estimands per dataset, bias, and more datasets chosen by a
  stated rule.

**Design.**
- Part A, the reread: below.
- Part B, no compute: Murray's a0 per real dataset, and where it falls under the default's induced prior on the
  spread of the log mean.
- Part C: arms the default, Murray's calibration fixed per dataset (k = 3/a0) and a drawn prior centered at
  Murray's value; the 11 real datasets with the landed study's two truths each (as fitted, and BART-fitted), 20
  replications; each dataset's estimand as its own literature reports it. Real outcomes by Murray's 10-fold
  log loss and WAIC; there, a cross-validated arm over k in {0.5, 1, 2, 4}, by a loop, since `xbart` does not
  take the family.
- If the maintainer adds datasets, the added ones run Part C at the same size.

**Readings.** Per dataset and estimand: RMSE, bias, coverage of 90 percent intervals and width, with intervals,
on its own spread, never pooled across datasets first; real outcomes as Murray's Table 3, log loss and WAIC per
dataset. No bar.

**What results mean.** The default's RMSE, bias and coverage per dataset alongside Murray's: where Murray's
calibration does as well or better on most datasets and never much worse, the recommendation is Murray's
calibration as the default; otherwise the default stays, with the datasets where it loses named.

**Stored results.** The landed study's scores hold, for the default at its default length, 5 replications per
real dataset (10 for crimtab), coverage of 50, 80, 90 and 95 percent intervals, interval score, width, RMSE
and log density by estimand, for each truth and its widened, subsampled and fold cells. Its full draws, with
truths, are archived on the laptop, so bias reads at no compute. A first pass: the default's 90 percent coverage
of fitted means runs from 0.76 (Traffic, BART-fitted truth) to 1.00 (Insurance), and of rate ratios from 0.79 to
1.00. The reread gives every dataset's own table with bias added; with 5 replications its intervals are wide,
so it reports and does not settle.

**Report outline.** 1 page question, Murray's calibration and Part B's figure; 2 pages per-dataset tables
(RMSE, bias, coverage, width) by estimand; 1 page real outcomes in Murray's form; 1 page meaning and limits.

**The study's own choices.** The drawn prior centered at Murray (alternative: Murray's fixed arm alone); the CV
grid (alternative: Murray's evaluation alone, no CV arm); 20 replications (alternative: 50, about 2.5 times
Part C's cost); the two truths (alternative: as fitted only).

## A3 with A4. stan4bart's initial fit and the starting sigma on every forest's basis

### A3. The initial fit

**Question.** A stan4bart user gets the forest's range from an initial linear mixed model that holds the forest's
columns as linear terms. Does that, against splines on those columns, change the treatment effects, group effects
and intervals they report, across designs and on real multilevel data?

**Published reports.**
- Dorie, Perrett, Hill and Goodrich (2022), sections 6 and 7: IHDP covariates, grouping by site (8) or by
  mother's age (26), varying intercepts or intercepts and slopes, response surfaces A, B and C as Hill (2011),
  intraclass correlation 0.2, 0.333 and 0.5, effect sizes 0 to 0.8; SATT, group ATT and PEHE, RMSE and interval
  length divided by the outcome's sd, 95 percent coverage. The replication count is not stated in the text read.
- Hill (2011): 1000 replications of IHDP surfaces.

**Design.** Arms: linear (shipped) and three-degree-of-freedom splines, the landed study's two. Cells: 12 of the
paper's IHDP settings (two groupings, two random-effect forms, surfaces A, B and C at ICC 0.333), 100
replications; the landed study's half-curved family, 100 seeds; radon (919 homes, 85 counties) and High School
and Beyond (7185 students, 160 schools) with forest-fitted truths, 100 and 50 replications; both datasets' real
outcomes by 10-fold cross-validation.

**Readings.** The paper's: SATT and group ATT RMSE over the outcome sd, PEHE, coverage and interval length; on
radon and HSB, county or school means and the fixed effects (log uranium; sector, mean SES), RMSE, bias and
coverage; real outcomes by held-out log density. No bar.

**What results mean.** Splines better on some cell by more than the landed study could resolve, and no worse on
the others: the recommendation is splines, or a documented knob; otherwise linear stands, with the cells named.

**Stored results.** The landed study's 240 in-between fits and its earlier stress fits are not in their summary
form; its tables remain. No reread beyond them.

### A4. The starting sigma on every forest's basis

**Question.** A user writing a two-forest treatment model gets a starting sigma from a regression on x and the
treatment basis. Does that, against the regression on x alone and against the sample sd, change the effect
estimates, their intervals and the noise sd's coverage?

**Published reports.**
- CGM (2010), section 2.2.4: the sigma prior is set from "a rough data-based overestimate": the naive sample sd
  of y, or the residual sd of a linear regression; q = 0.90 by default.
- The bcf package: the residual sd of `lm(y ~ z + x)`, nu 3, q 0.9.
- Hahn, Murray and Carvalho (2020), section 6.1: 8 designs (n 250 or 500, homogeneous or heterogeneous effect,
  linear or nonlinear prognostic function), 200 replications, ATE and CATE RMSE, coverage and interval length
  (Tables 2 and 3).

**Design.** Arms: shipped (x and the basis, bcf's choice), x alone (the old default, by an explicit `sigest`) and
the sample sd (CGM's naive choice), the last at 50 replications. Cells: BCF's 8 designs at 200 replications;
IHDP surface B at 100.

**Readings.** BCF's Tables 2 and 3: ATE and CATE RMSE, coverage and length of their intervals; beside them,
coverage of the noise sd's interval and the share of replications where the starting estimate is at least the
true sd (CGM's intent of an overestimate). No bar.

**What results mean.** The effect readings the same across arms (as the landed 20-seed check found) and the noise
sd's coverage worse under shipped at small n: the report says so, and the recommendation is to keep bcf's choice
and document the cost; a difference in effect intervals is a recommendation to revisit.

**Stored results.** The landed 20-seed check is summarized in its design note only. No reread.

**Report outline (A3 and A4 together).** 1 page questions and published designs; 1.5 pages A3 in the paper's
tables; 1 page A4 in BCF's tables; 1 page radon and HSB; half a page limits.

**The study's own choices.** A3: 12 of the paper's settings at one ICC (alternative: all three ICCs, about 3
times the IHDP cost); 100 replications (alternative: Hill's 1000); forest-fitted truths for radon and HSB
(alternative: real outcomes only). A4: the overestimate share (alternative: noise sd coverage alone); the naive
arm at 50 replications (alternative: 200).

## A6. The monotone defaults

**Question.** A user constraining f to be monotone gets the shipped tree prior. Do 95 percent intervals for f
cover near nominal, and is the fit's tree size what an unconstrained fit gives, as mBART intends?

**Published reports.**
- Chipman, George, McCulloch and Shively (2022), mBART, section 4.3: alpha .25 and beta .8, set to make up for an
  approximation in their sampler, aiming at "tree sizes comparable to those obtained in unconstrained BART".
- mBART section 5.3: f = x1 x2^2 + x3 x4^4 + x5 on five uniform predictors, sigma in {0.2, 0.5, 0.7, 1, 2}, 200
  datasets of 500 training and 1000 test rows; RMSE ratio boxplots (Figure 8) and Table 1's width and coverage;
  mBART covered 90.1 to 97.9 percent, read as "good reliability", with no bar.

**Design.** Arms: unconstrained, the shipped prior with the constraint, and mBART's prior with the constraint.
- C1, the reread: below.
- C2, the partial rerun, at 75 trees and 4 chains: the stored cells the reread leaves open (the three low-noise
  cells at n 2000), 50 replications, and mBART's section 5.3 at its five noise levels, 50 replications.

**Readings.** As mBART's Table 1 and Figure 8: coverage and width of 95 percent intervals for f, RMSE ratios
against the unconstrained fit, and splits per tree against the unconstrained fit; with intervals. No bar.

**What results mean.** The shipped prior's tree size nearer the unconstrained fit's than mBART's prior is, and its
coverage between the two: the shipped prior meets mBART's own aim under an exact sampler, and stands.

**Stored results.** The tree-prior study's baseline holds 24 cells (three designs, two noise levels, n 200 and
2000, 50 and 200 trees), 8 replications, one chain of 500 after 500. Reread: in all 24 cells coverage orders
unconstrained, then shipped, then mBART's prior; the shipped prior's splits per tree are 1.27 to 1.32 against the
unconstrained fit's 1.40 to 1.45 and mBART's prior's 0.29 to 0.68. So by mBART's own criterion the shipped prior
is the right one, at both tree counts that bracket the default. Not settled: at 50 trees, low noise and n 2000
every arm covers below nominal (the unconstrained fit 0.84 to 0.87), a property of the forest and not of the
constraint, and the default's 4 chains are not in the baseline. The settle rule (the project's own): the reread
settles a cell when the ordering of the arms is the same at 50 and 200 trees; the alternative is always running
the rerun.

**Report outline.** 1 page question and mBART's calibration; 1.5 pages the reread's table and tree-size figure;
1.5 pages C2 in mBART's Table 1 and Figure 8 form; 1 page meaning and limits.

**The study's own choices.** The settle rule; 50 replications on mBART's design (alternative: its 200, about 4
times C2's cost); the open cells chosen by the reread (alternative: all 12 cells at 75 trees).

## Sizing

Costs are core-hours on the box; precision is the half-width of a 95 percent interval at full size, estimated
from the stored spreads; resolution targets are the project's own (the landed studies' resolution).

| study | component | question it serves | size (looks) | expected | at most | precision at full size | resolution target | stop |
|---|---|---|---|---|---|---|---|---|
| A5 | C1 | claim on the landed designs | 5 cells x 2 arms x 25 (13, 25) | 0.2 | 0.3 | coverage 0.01, RMSE 1% | 0.02, 1.2% | at 13 if every interval is inside target |
| A5 | C2 | claim on Tan's data | 12 cells x 2 x 25 (13, 25) | 0.4 | 0.6 | coverage 0.01 | 0.02 | as C1 |
| A5 | C3 | mixing at Tan's length | 3 cells x 2 x 10 (5, 10) | 2.0 | 3.0 | R-hat 0.005 | none: reported | full size |
| A5 | C4 | 200 trees | 5 x 2 x 25 (13, 25) | 0.6 | 0.9 | as C1 | as C1 | as C1 |
| A5 | C5 | recovery after a response swap | 2 x 2 x 25 (13, 25) | 0.2 | 0.3 | 10% of sweeps | 20% | as C1 |
| A2 | A, B | prior; reread | none | 0 | 0 | | | |
| A2 | C | default against BART-cv | 22 datasets x 20 splits (10, 20) | 2.7 | 4.0 | median relative loss 0.005 | 0.01 | at 10 splits if inside target |
| A2 | D | calibration on truths | 18 cells x 3 x 50 (25, 50) | 3.5 | 5.3 | coverage 0.01 | 0.02 | at 25 if inside target |
| A1 | A, B | reread; Murray's a0 | none | 0.1 | 0.1 | | | |
| A1 | C, truths | per-dataset estimands | 22 cells x 3 x 20 (10, 20) | 2.0 | 3.0 | coverage 0.03 | 0.05 | at 10 if inside target |
| A1 | C, folds | real outcomes | 11 x 10 folds x 4 arms | 1.7 | 2.5 | log loss 1% | none: reported | full size |
| A3 | IHDP | the paper's settings | 12 x 2 x 100 (50, 100) | 2.0 | 3.0 | coverage 0.02 | 0.03 | at 50 if inside target |
| A3 | others | half-curved, radon, HSB, folds | as designed | 1.2 | 1.8 | split error 2% | 3% | as IHDP |
| A4 | BCF | effects and sigma | 8 x 2 x 200, naive x 50 (100, 200) | 1.4 | 2.1 | coverage 0.02 | 0.03 | at 100 if inside target |
| A4 | IHDP B | effects on real covariates | 3 x 100 (50, 100) | 0.3 | 0.5 | as BCF | as BCF | as BCF |
| A6 | C1 | reread | none | 0 | 0 | | | |
| A6 | C2 | open cells, mBART's design | 8 cells x 3 x 50 (25, 50) | 5.0 | 7.5 | coverage 0.01 | 0.02 | at 25 if inside target |

Expected 23.3 core-hours (A5 3.4, A2 6.2, A1 3.8, A3 3.2, A4 1.7, A6 5.0), at most 34.9; added datasets for A1
cost about 0.5 core-hours each. Timings used: a default gaussian fit 1.8 to 2.7 s at n 500 to 2000, probit 1.6 s
at n 2000, logistic about 2, 4 and 8 s at n 500, 2000 and 4000, negative binomial a median of 4 s (2 to 25 s by
dataset), stan4bart on radon 2.6 s, a two-forest treatment fit 1.2 s at n 500, a constrained monotone sweep 4 to
5 ms at 50 trees.

## Calls for the maintainer

1. A1's added datasets. (a) Independent agents choose about 12 count datasets by a stated rule, as dec-B635 does
   for pricing; (b) the maintainer names them; (c) the 11 alone. Recommend (a): the referees asked for it, and a
   rule is what published comparisons state.
2. The run order inside the tier. Recommend the order above, with A2's Part D after pricing's ruling and A6's
   rerun only on the cells the reread leaves open.
3. stan4bart 0.0-13 as an arm in A3. It does not load beside today's libraries; a private build against
   dbarts 0.9-34 on the box might. Recommend one attempt, the arm added if it loads and omitted, with that said,
   if not.
4. Where a study recommends a change of default or of a documented claim (A5's news entry, A2's k, A1's
   calibration), that change is the maintainer's; each report recommends one.
5. The budget: about 23 core-hours expected, about 35 at most; about 29 and 44 with 12 added datasets for A1.
   Recommend approval per study, in order.

## Papers

Read in full: Chipman, George and McCulloch (1998, 2010); Chipman, George, McCulloch and Shively (2022); Tan,
Ronen, Saarinen and Yu (2024); Murray (2021, the arXiv version); Pratola (2016); He and Hahn (2023); Hahn, Murray
and Carvalho (2020); Hill (2011); Hill, Linero and Murray (2020); Vehtari et al. (2021); Dorie, Perrett, Hill and
Goodrich (2022, the full text from Europe PMC; the publisher's PDF was refused). Not reached: mBART's used-car
data (section 5.4), whose availability is unverified; the stan4bart paper's replication count, not in the text.

## Verification

- Run: 0.9-34's tree-move mixture and binary k prior on the released package; the induced logit prior under
  the default, fixed k = 2 and Murray's k; the reread of the k study's real splits and of the monotone baseline
  quoted above; the first pass of the negative binomial reread; the timings in [Sizing](#sizing).
- Read: the bcf package's starting sigma; stan4bart 0.0-13's initial fit, in its source.

## Not verified

- Box timings: every cost is from the laptop under load, to be rechecked in each smoke.
- That PMLB's copies of Tan's four datasets match Tan's, and that Tan's synthetic designs are copied exactly.
- The precision column: estimated from stored spreads, not from a pilot.
- That the archived negative binomial draws read back whole.
