# binary-k-mixing-need: does slow k mixing harm the estimands users report, and where?

Status: PLANNED 2026-10-17, PRE-REGISTERED 2026-10-29, amended the same day after review. Runs after the
merge to main and before 1.0-0 (dec-B462), behind the work that stops the merge. The study dec-B462 rules,
revising dec-B456 for it; its bar is dec-B455's; its use tiers count a logistic cell as the same cell
under probit (dec-B463). Revised from its blind critique (scratch/kmn/crit/crit-report.md, SOUND AFTER
CHANGES) and the orchestrator's calls on it ([Calls made in planning](#calls-made-in-planning)). Stage 0,
stage 1 and the smoke ran 2026-10-29; the rule, bands and bars are fixed in
[Pre-registration](#pre-registration). No study fit yet. One slot waits on the maintainer's colleagues
([Rulings and pending slots](#rulings-and-pending-slots)).

agent: one opus implementer writes the harness, runs it on the bench box and writes the finding; one
review of the harness and smoke (the gate-only tier of the plan README); one review of the tranche 1
reading before tranche 2 runs (opus). The blind critique ran and is answered here.
rng: NEUTRAL. No package code changes: a benchmark harness, two private builds of the collapsed move's
branch and stan4bart installed against each.
window: after the merge to main, before 1.0-0 (dec-B462). Its result decides whether the collapsed move,
or any faster k move, is priced at all.
budget: harness about 1,900 lines of R, much of it reused from the binary interval study and the SBC
harness. Compute about 141 core-hours expected (about 40 hours of the box), 167 at most, stop at 211
([Cost](#cost)).

## Summary for the maintainer

You ruled that k must move for binary fits (dec-B420) and asked the next question: is any movement
enough, or does k have to move fast and reach its stationary distribution for the estimands people
report to come out right? And if fast moves are needed, where, and does that matter for how dbarts is
mostly used?

What is known: at the default length the shipped sampler's k is unconverged in nearly every binary fit,
under both links; the collapsed move mixes k about a hundred times faster. But the pilot of the binary
interval study found no gain from the move in treatment effects, false declarations or interior
probabilities, and a corrected re-read of its real-outcome reading found both samplers reproducing the
long-run posterior about equally well where the long run could be trusted (0.895 against 0.900 under
logistic, against a nominal 0.90). The only places a slow k has been seen to matter are probabilities very
near 0 or 1, the reported k itself, and regimes nobody has tested.

This study tests them. It asks, regime by regime and estimand by estimand, whether the shipped sampler
at the default length gets the location and uncertainty of what users report wrong by at least 2 points
of interval coverage, and whether the error is k's. Two measures decide:

- **Calibration against a known truth.** The truth is drawn from the model's own prior on real
  covariates, so a sampler that reaches the posterior covers at nominal by construction, and any gap is
  the sampler's (simulation-based calibration, SBC). A companion fit with k fixed at its true value
  shows what the trees alone achieve, so whether a gap is k's can be read off.
- **Agreement with a long-run posterior** that both samplers reach and agree on: the share of the long
  run's draws inside a default-length interval, against the 0.90 a 90 percent interval should hold.

Fits with k fixed at quantiles of its posterior bound what a k that never moves can do; pooled, they
mimic a run whose chains each sit still at a k drawn from its posterior. The ends of the long runs give
each kernel started at stationarity, which separates harm from a far start (fixable cheaply, for example
by moving k only during burn-in) from harm from slow mixing.

Regimes: near separation that is not deterministic, rare outcomes (1 to 5 percent, weak signal) at 500
and 10,000 rows, small samples with strong signal, many predictors, a moderate control, treatment effects
at moderate and rare rates, the propensity stage with propensities used as weights, your
multiple-treatments design, real outcomes, and BART inside stan4bart with a varying intercept and inside a
dbartsSampler loop. Estimands: the probability (interior and near 0 or 1), log risk and odds ratios at
rare rates, average effects, individual and subgroup effects, the multiple treatments' contrasts, false
declarations and bands, a propensity-weighted effect, inverse-propensity weights, and the stan4bart group
standard deviation.

Each cell ends as "problem shown", "no problem" or "unresolved", and cells are weighed by how common their
use is: effects at moderate prevalence, propensity weighting, your colleagues' case and embedded use
first; rare outcomes, many predictors and small samples next; near-separated prediction last; a logistic
cell counts as the same cell under probit (dec-B463), since shipped probit has a rescaling step for k and
logistic has none, so logistic is where a problem is likeliest. A problem shown where
use is common opens the pricing of moves; a problem shown only in rarer use leads to a scoped or opt-in
option; no problem is written down with triggers for reopening, and the move stays parked on its branch.
If no problem is shown but some common-use cells stay unresolved, the move stays parked and those cells
come to you with what resolving each would cost.

Cost: about 141 core-hours expected, about 40 hours of the box, in two tranches; the second runs only
where the first leaves a cell outside near-separated prediction unresolved or shows a problem. The stop is
at 211.

## Goal

After this study each regime and estimand class of [Regimes](#regimes) and [Estimands](#estimands) carries
one reading per link, fixed by the rule of [What counts as a problem](#what-counts-as-a-problem): problem
shown, no problem (bounded or measured), or unresolved, with the effective draws per second of every
estimand under both kernels and the share of any harm due to the start. That map, weighed by use, says
whether faster k moves are worth pricing and where.

## What is known

- Binary fits draw k under chi(1.5, 2) by default (dec-B420,
  [binary-k-fixed-or-drawn.md](../design/binary-k-fixed-or-drawn.md)). The best fixed k lost on real data
  in both directions, so k must adapt; what was not settled is whether it must mix. Under logistic k is
  unconverged at 64,000 draws while the predicted probabilities and log score are converged; the reported
  k and intervals very near 0 or 1 depend on chain length.
- How k enters: each leaf's prior scale goes as 1/k, so k sets the overall scale of the latent fit.
  Interior probabilities are pinned by the data and move little with k; near separation only the prior
  holds the scale back, its posterior is wide and slow, and the tail probabilities follow it
  ([k-mixing-pg-families.md](../design/k-mixing-pg-families.md)). A functional's effective draws depend
  on how much of it lies along that slow direction, not on k's own autocorrelation.
- The collapsed move (leaves times a factor, k divided by it, the factor drawn with the latents
  integrated out, one sweep in four) is built on the branch wt/collapsed-move-plan, tip 207f740e
  ([k-tail-samplers.md](../design/k-tail-samplers.md) for the survey). Its cost on the box, from the
  binary interval study's timing pass: 1.006 to 1.30 times a shipped iteration, mostly 1.03 to 1.15.
  Even with it, k's split-Rhat at the default length has a median of 1.1 to 1.25 in the pilot's cells and
  reaches 1.6 on the effect designs (the MCMC panel's reading).
- The binary interval study ([binary-interval-study.md](binary-interval-study.md)) stopped after its
  pilot (dec-B462). Its pilot: no gain from the move in effects (T +0.004 probit, +0.009 logistic),
  false declarations or the multiple-treatments cells; its apparent harm on probabilities sat entirely
  in near-deterministic cells, where the truth, a posterior mean, lay outside the posterior's own
  interval. The corrected real-outcome re-read (f30d5c82 on its branch) put agreement with the long run
  on interior rows at 0.872 against 0.881 (logistic) and 0.878 against 0.876 (probit), and on accepted
  references at 0.895 against 0.900 and 0.898 against 0.898, against a nominal 0.90. The earlier 0.75
  against 0.86 was saturated rows and single-precision rounding; no arm here rests on it. Only 18 of its
  52 references passed its acceptance, which required k itself to converge; most failed on k in the
  near-separated and large datasets, exactly where the question is open.
- The panel convened on the pilot (MCMC, applied, design, maintenance) named the untested ground: ratio
  estimands at rare rates, propensity models near positivity and propensity weighting, rare outcomes at
  large n, simultaneous bands outside near-deterministic cells, the multiple-treatments design, and BART
  inside a larger sampler. It proposed the fixed-k quantile arms, the variance of each estimand explained
  by k, chains started apart, SBC truths, a typical draw instead of the posterior mean as truth, and
  reading by regime.
- The default link: bart() and bartCause take a 0/1 outcome as probit unless asked otherwise, and
  stan4bart fits a binary outcome only under probit. The samplers differ more than the models: shipped
  probit carries a rescaling step for k (its 99th percentile of k's autocorrelation time about 900 sweeps,
  from about 23,000 without it) and logistic carries none (about 22,600), so logistic is where a problem
  is likeliest, and the maintainer ruled that a logistic cell counts as the same cell under probit
  (dec-B463).
- stan4bart's BART part is a dbarts probit sampler taking one sweep per outer iteration, drawing k under
  the dbarts default unless its BART arguments fix it, beside a parametric part (fixed and varying
  effects) sampled by WALNUTS with an exact move on each group scale. It has no chain continuation; it
  offers a thinning setting and a per-iteration callback. The probit half of the collapsed move reaches
  every such fit that draws k; the logistic half reaches none. A probit 0/1 weight there is the
  active-row mask, so all-zero weights sample the prior; its default missing-value handling drops rows.
- The SBC harness [benchmarks/R/sbc.R](../../benchmarks/R/sbc.R) draws a truth from the sampler's own
  prior, k included, through [`sbcReplication`](../../benchmarks/R/sbc.R) and
  [`sbcKDraw`](../../benchmarks/R/sbc.R), with self-checks
  ([`sbcCheckKHandBack`](../../benchmarks/R/sbc.R), [`sbcCheckFitConsistency`](../../benchmarks/R/sbc.R)).
  It covers probit with a drawn k on simulated covariates today.

## The question, made operational

A slow k can harm a reported estimand in three ways, and the arms are chosen to separate them:

1. **Stuck but stationary**: k sits in its posterior but does not travel through it in one run, so the
   run reports the estimand's distribution at one k, or at a few, rather than averaged over k. Bounded by
   the fixed-k arms and their pooled stuck runs.
2. **Transient**: the chain starts outside the bulk and has not arrived after the default burn-in, so
   the run reports the wrong location. Measured by each kernel's default against the same kernel started
   at stationarity (the reference tails); k's share of it by comparing the shipped kernel's transient with
   the move's from the same start, since the trees also start from stumps.
3. **Slow at stationarity**: the run is in the right place but has few effective draws of what moves
   with k. Measured by the two kernels' stationary tails against each other, and by effective draws per
   second.

"Any movement is enough" for an estimand means the three together cost less than the bar. "Fast moves
are needed" means they cost at least the bar and a sampler that mixes k faster removes it.

## What counts as a problem

Fixed before any result. Each reading is per link, per regime of [Regimes](#regimes), per estimand class
of [Estimands](#estimands).

### Measures

SBC cells:
- **G**, the shipped default's coverage gap: the class's 90 percent pointwise coverage of the truth,
  pooled over all replications first (SBC is nominal only on average over them), then its distance from
  0.90 on either side, corrected for noise as dec-B455 defines it, with the noise variance taken from the
  bootstrap over replications rather than from binomial counts (rows within a replication share one
  truth). 95 percent read beside.
- **D**, the signed paired difference in pooled coverage, shipped default minus the oracle (k fixed at
  its true value, nominal by construction once its trees mix; [Arms](#arms)). Signed, so that errors on
  opposite sides of nominal cannot cancel as a difference of distances would.

Reference cells (every regime's reference replications, and D1, D2 and B):
- **Total shortfall S**: 0.90 minus the shipped default's reference mass, the share of the reference's
  draws inside its 90 percent interval, averaged over the class's quantities and over replications. The
  comparator is nominal: the interval between the order statistics of N independent draws holds, in
  expectation, exactly its nominal share of the posterior, and a stationary chain departs from that only
  by a term of order one over its effective draws; the stationary tails show how large that is in each
  cell.
- **Breakdown of S** (pricing reads it; the problem reading does not need it): the shipped transient
  (shipped tails minus shipped default), the move's transient (the same for the move), k's share of the
  transient (the first minus the second), and k's share at stationarity (move tails minus shipped tails).
- **Stuck-run shortfall**: 0.90 minus the reference mass of a pooled stuck run (below), and, for the
  ruled fixed-k arms at q10, q50 and q90, 0.90 minus each one's mass.
- Every reading of a kernel is scored against the other kernel's half of the reference: the shipped
  default and the shipped tails against the move half, the move's against the shipped half, so a tail is
  never scored against its own chain and each transient is paired on one half. The pooled reference less
  the scored chains is read beside, as a check.

Secondary, reported, not deciding: the location error, |default posterior mean - reference mean| over the
reference standard deviation; the width ratio to the reference interval; the probability-scale interval
score (Gneiting and Raftery 2007); under SBC the dispersion ratio, the mean over replications of (truth -
posterior mean)^2 / posterior variance, which is exactly 1 for an exact sampler; and, for designed truths
near 0 or 1, coverage with the interval widened by 0.001 on the probability scale.

### Bars and where they come from

| quantity | bar | source |
|---|---|---|
| G, D, S, stuck-run shortfall, transients | 0.02 at 90 percent | The house's noticeable coverage difference: September's binary hyperprior rule, reused by the k study ([The bars](../design/binary-k-fixed-or-drawn.md#the-bars)) and ruled for calibration in dec-B455. Its meaning: for a near-normal posterior a 90 percent interval that covers 0.88 is 5.5 percent too narrow, and its miss rate is a fifth higher. |
| the same, read in reference mass | 0.02 | An interval holding a share m of the posterior covers a truth drawn from that posterior with probability m, so mass and coverage share one scale. |
| location error | 0.35 posterior sd | Derived: for a normal posterior, a shift of 0.35 sd costs a 90 percent interval 0.02 of coverage. |
| width ratio | 0.945 or 1.06 | Derived the same way: the scaling that costs 0.02. |
| tail rows | probability outside [0.01, 0.99] | The house's interior split (k study, binary interval study), applied to a data-defined label ([Truth construction](#truth-construction)). |

One bar for every class, because the bar states what a user would notice, not what the counts can
resolve; where the counts cannot resolve it, the reading is "unresolved", not a looser bar. The binary
interval study's false-declaration bars (M1 0.02, M2 0.05, M3 0.01) are not used to decide: they were gain
bars between two samplers with no comparator here. D1's false declarations are read as reference mass of
the contrasts that decide them and reported as rates.

### Bounds

Every bound is the estimate plus or minus a quantile times its standard error, the standard error from a
bootstrap of 4,000 resamples over replications with all cells of a replication resampled together (4,000
resamples cannot place a tail quantile this far out, but they estimate a standard error well). The
quantile is normal for SBC readings and Student's t with one degree of freedom fewer than the cell's
passing references for reference readings. Levels:
- **problem**: one-sided 0.000156, z = 3.60: Bonferroni over about 80 cells a link, two deciding
  measures a cell and two looks (320 tests), so the chance of any false "problem shown" in a link is under
  0.05; t at the same level (about 4.4 at 24 references).
- **no problem**: one-sided 0.025 a look, z = 1.96 (t 2.07 at 24 references, 2.20 at 12).

### The reading of a cell

- **Problem shown**, either of:
  - (SBC) G at least 0.02 with its lower bound above 0, and D at least 0.02 in size with its bound
    excluding 0 on the same side. A gap with D inside the bar is the trees' or the model's and goes to the
    diagnosis, not to k.
  - (references) S at least 0.02 with its lower bound above 0, and k's by either test: the move's default
    from the same start is within the bar of 0.90 while the shipped default is not, or the stuck-run
    shortfall is itself at least 0.02 with its lower bound above 0.
- **No problem, bounded** (references): the stuck-run shortfall and k's share of the transient both have
  their upper bound below 0.02. A run whose chains each sit at a stationary k, anywhere in k's posterior as
  seven quantiles represent it, cannot cost the bar, and the default start adds less than the bar. The
  first half does not depend on how the shipped kernel moves k; the second does, and is measured.
- **No problem, measured**, all available of: (SBC) G's upper bound below 0.02, or D's two-sided bounds
  inside (-0.02, 0.02) (then any gap is not k's); (references) S's upper bound below 0.02.
- **Unresolved**: anything else, with the reason (counts, failed references, or both).
- A cell's reference readings count only when at least 80 percent of its attempted references pass
  acceptance, and only if the reading does not change when each failing replication is scored against its
  move half alone. Otherwise its reference readings are unresolved, and the reason is stated.

The problem readings use the shipped sampler at the default length, since that is what users run;
equal-time comparisons belong to pricing.

### Reading by regime and weighing by use

The result is a map: rows are regimes, columns estimand classes, one map per link, each cell one of the
readings, with the breakdown of any shortfall. A regime needs efficient k moves when any estimand cell in
it shows a problem that is k's.

Cells are weighed by use in three tiers, lexicographically (a problem in a more common tier is never
outweighed by its absence elsewhere). The tiers follow the applied panel's ordering of what users report
and where they work, the maintainer's "more important estimands", the consumer map (bartCause, stan4bart
and treatSens fit causal effects; stan4bart and dbartsSampler loops are dbarts's distinguishing use). A
cell's tier is the rarer of its regime's and its estimand's. The link does not enter: a logistic cell
counts as the same cell under probit (dec-B463), because shipped probit has a rescaling step for k and
logistic has none, so logistic is where a problem is likeliest; a fix a logistic-only finding leads to
still touches only logistic fits.

| tier | regimes | estimands |
|---|---|---|
| 1, the bulk of use | moderate control; effects at moderate prevalence; propensity weighting near positivity; the multiple-treatments design and the colleagues' case; stan4bart and the dbartsSampler loop | average effects (overall, treated), subgroup averages, interior probabilities, the multiple treatments' contrasts, false declarations and bands, the weighted effect and inverse weights, the group standard deviation |
| 2, common in health data | rare outcomes at small and large n; many predictors; small samples with strong signal | individual effects, ratio estimands at rare rates |
| 3, niche | near separation | probabilities near 0 or 1, other simultaneous bands, k |

Each real dataset of B takes the tier of the regime its separation score and outcome rate place it in.
What follows from each tier is under [What the result leads to](#what-the-result-leads-to).

## Regimes

Real covariates throughout; outcomes, and for S9 the treatment, are drawn. Selection uses the data alone
([Truth construction](#truth-construction) shows why coverage stays nominal). Links: both in every SBC
prediction regime, S7, S8, E3, D1 and B; probit only in S9, D2, E1 and E2 (below).

| id | regime | covariates, n | selection on the data | reps (tranche 1, full) | references (1, full) |
|---|---|---|---|---|---|
| S1 | near separation, not deterministic | tictactoe, banknote, wdbc, ionosphere; 500 or all rows | separation score in the near-separated band | 100, 400 | 12, 24 |
| S2 | rare outcome, small n | pima, german, creditapproval; 500 | outcome rate in [0.01, 0.05] or [0.95, 0.99] and separation score at or below the moderate band | 100, 400 | 12, 24 |
| S3 | rare outcome, large n | adult, 10,000 subsampled | as S2 | 40, 80 | 6, 8 |
| S4 | small n, strong signal | kyphosis, birthwt, haberman, cleveland, infert; 81 to 200 | strong band, below near-separated | 100, 400 | 12, 24 |
| S5 | many predictors | ACIC 2016's 58 covariates with 192 permuted copies, p = 250; 1,000 | none | 100, 200 | 12, 20 |
| S6 | moderate control | pima, german; 500 to 768 | moderate band, rate in [0.2, 0.8] | 100, 400 | 12, 24 |
| S7 | effects, moderate | IHDP, 747, its real treatment | none | 100, 300 | 12, 24 |
| S8 | effects, rare outcome | ACIC 2016, 1,000, setting 31's treatment model | as S2 | 100, 200 | 12, 20 |
| S9 | propensity near positivity, used as weights (probit) | IHDP covariates, 747; treatment drawn | weak overlap by an auxiliary propensity | 100, 300 | 12, 24 |
| E1 | stan4bart, varying intercept (probit) | Twins, 1,000 complete pairs subsampled, birth state as group (about 50) | rate in [0.1, 0.9] | 100, 200 | 12, 20 |
| E2 | stan4bart, few groups, rare (probit) | as E1, states pooled into 8 regions | rate in [0.02, 0.08], separation as S2 | 100, 200 | 12, 20 |
| E3 | dbartsSampler loop: 30 percent of outcomes imputed every sweep | pima, german; 500 | moderate band | 100, 200 | 12, 20 |
| D1 | multiple treatments (the binary interval study's A6) | ACIC 2016, 2,000, four arms; its four cells | designed truth | reference replications only; its pilot's 10 a cell as a secondary reading | 4 a cell |
| D1x | the colleagues' case, if shared | theirs | as given | as D1 | 4 to 6 |
| D2 | propensity-weighted effect, weak overlap (probit) | ACIC 2016 setting 1, 1,000; base rate 0.3 and 0.05 | designed truth | reference replications only | 6 a cell |
| B | real outcomes | the 22 real datasets and adult (4,000 rows) of the k study; lindner, ACTG 175, lalonde | none | 6 seeds | 1 a problem |

Notes on the regimes:
- **Separation score**: the cross-validated AUC (5 folds) of a classification tree (rpart, a recommended
  package) of the outcome on the covariates; a function of the data alone. Its bands are set in stage 1,
  before any arm runs, from the k study's simulated cells whose truth-based sign agreement is about 76
  percent (moderate) and 90 percent (near-separated), the definitions the k study and the collapsed
  plan's large cells used: the moderate band and the near-separated band are the central 80 percent of
  those cells' scores, the strong band lies between them. A score of 0.995 or more is deterministic and
  rejected: nobody reports intervals there and the data cannot tell 1e-4 from 1e-8 (the binary interval
  study's near-deterministic finding). Results of every regime are also reported by separation score.
- **Rare outcomes, weak signal**: with the default offset of 0, a rate of 1 to 5 percent needs a large
  negative shift shared by all rows, and draws that reach it mostly carry a large spread across rows: a
  few near-certain positives and the rest near 0, which is quasi-separated small-k data, not the
  weak-signal rare outcomes of health data. S2, S3, S8 and E2 therefore also require the separation score
  at or below the moderate band, computed only on draws that pass the rate test, so it stays cheap. Both
  tails are accepted: with offset 0 the prior is symmetric under turning f to -f and y to 1 - y, under both
  links. Acceptance may be well under 1 percent; what is gated in stage 1 is the cost per accepted
  replication, which may not exceed the replication's fit cost (prior draws are cheap; the tree score runs
  only on rate passes). A regime over that gate runs at the replications it affords, and says so.
- **Large n** is 10,000 rows: the collapsed plan's large-cell size and the low end of the electronic
  health records the applied panel names. 20,000 or 50,000 rows double or quintuple S3's cost.
- **S5's truth is dense**: the BART prior splits uniformly over all 250 covariates, so the truth depends
  on the permuted copies as much as on the real ones. S5 tests posteriors with many predictors, not sparse
  signal, which would need a sparsity prior both in the truth and in the fit.
- **S7 to S9** follow bartCause's default response model: the outcome on the covariates, the treatment
  and the estimated propensity, 10 chains. In S7 and S8 the propensity is the shipped build's 10-chain
  probit fit, made once per design before the truth is drawn, so every arm faces the same design. S9
  draws the treatment from a probit BART prior propensity and fits the propensity with each arm; because
  the truth cannot be drawn on an arm's own propensity, its outcome model leaves the propensity out (a
  bartCause option), which keeps the joint posterior exact (see [Truth construction](#truth-construction)).
  D2 keeps bartCause's default, the propensity as a covariate, with each arm's own propensity, at the
  price of a designed truth.
- **S9 and D2 are probit only**: bartCause fits its propensity under probit whatever the outcome's link,
  and S9's prior propensity is probit, so the propensity stage users meet is probit even in a logistic
  outcome fit. A logistic half would cost about 5 core-hours in tranche 1 and 6 in tranche 2.
- **Weak overlap** (S9): at least 10 percent of units with an auxiliary logistic regression propensity
  outside [0.05, 0.95], a data-only rule; 0.05 is where bartCause's default bounds (0.025) start to bind
  within a factor of two.
- **E2** pools groups to 8 because stan4bart's own mixing note found the group scale weakest with few
  large groups, and a rare outcome because that is where the latent scale is least identified and k and
  the group scale can trade off. E1 and E2 draw their truth on the rows stan4bart keeps, the pairs complete
  on the covariates used.
- **E3** is the swapped-response use: a dbartsSampler run one sweep at a time, the 30 percent of outcomes
  missing at random redrawn from the current fit and set as the response between sweeps. Its target is
  the posterior given the observed rows, which plain fits on those rows also target, so the reference is
  made by plain fits.
- **D1 and D2 are reference cells.** A designed truth has no oracle and no nominal guarantee, so they are
  read only by the reference readings of [What counts as a problem](#what-counts-as-a-problem) (total
  shortfall, stuck runs, transient), on their reference replications, with every arm run on those. D1's
  pilot replications 1 to 10 (shipped and move, coverage, M1 to M3, bands) are a secondary reading
  against the designed truth, reported and not deciding; their stored draws are single-precision
  probabilities, which is harmless at D1's base rate of 0.3 except in its near-separated cell, where it is
  flagged. D2's designed-truth coverage on its reference replications is likewise secondary.
- **D1x**, the colleagues' multiple-treatments case: the D1 harness takes an external design and data; if
  shared before tranche 1 it runs as a tier 1 reference cell with 4 to 6 references (about 5 core-hours,
  within the stop), and if later, as an addition.
- **B** reuses the binary interval study's datasets, splits, pilot fits and time ratios. It has no truth;
  it reads the reference measures on real outcomes, where users are. adult at full size is dropped: its
  reference failed at 16 chains of 20,000 and costs about 1.7 core-hours a link to retry.

Not included, and why: deterministic outcomes (above); a multiple-treatments SBC (conditioning on which
contrasts are truly zero conditions on the truth, which SBC cannot score); the binary interval study's D3
cell, a rare effect at 4,802 units (cut on the critique's advice: S3 and S8 cover rare outcomes at large n
and in effects by SBC, which decides, where D3's designed truth could not); ACIC 2019 and Twins effect
cells (the pilot found nothing there and S5, S7, S8 and E1 cover their regimes).

## Estimands

Every interval is equal-tailed, of the per-draw quantity, at 90 and 95 percent, computed on the latent
scale where the quantity is a probability (equal-tailed intervals transform with a monotone map, so
coverage is the same and no rounding enters).

- **f**: the probability at every training row and at 100 held-out rows (S1 to S6, E1 to E3, B), split
  into interior and tail rows by a data-defined label.
- **Ratio estimands**: the marginal log risk ratio, log(mean p1 / mean p0), and log odds ratio of the
  marginal risks (S7, S8, D2); in the prediction regimes S2 and S3, the log risk ratio between the top and
  bottom quartiles of one covariate fixed per dataset in advance.
- **Average effects**: the sample average effect overall and on the treated, risk-difference scale (S7,
  S8, S9, D1, D2, B causal).
- **Individual effects**: each unit's p1 - p0; **subgroup averages**: four to ten subgroups per design
  defined by covariates in advance (for IHDP the birth weight quartiles and the mother's education levels).
- **Multiple treatments** (D1): the six pairwise contrasts per unit and their averages; whether each
  contrast is declared (its 95 percent interval excludes zero), read as the reference mass of the interval
  that decides it and reported as M1 to M3 rates with the binary interval study's definitions; the band
  over a unit's six contrasts and over all units', read as simultaneous reference mass.
- **Simultaneous bands**: the rank-based band of Besag et al. (1995) over held-out rows (f), over units
  (individual effects) and over D1's contrasts, built on the latent scale. Under SBC their coverage is read
  against the oracle at the same draw count, since a band built from order statistics depends on it.
- **Weighted effects**: bartCause's propensity-weighted effect on the treated, the sum of e_i times the
  unit's effect over the sum of e_i per draw (S9, D2). Because those weights are e_i, not 1/e_i, they barely
  feel the propensity's tails, so S9 adds the inverse-weighted (Hajek) estimates of the effect overall and
  on the treated, computed per draw from the propensity draws and the observed outcomes, with bartCause's
  default bounds on the propensity: a function of the posterior and the data, so still a posterior quantity
  under SBC. Also in S9 the propensity at each unit, the bounded inverse weights and the weights' effective
  sample size; in D2 bartCause's TMLE estimate and the largest inverse weight.
- **stan4bart** (E1, E2): the group standard deviation, each group's intercept, and the probability
  including group effects at each row.
- **k**: reported per arm, never read as an estimand: no user reports it.

## Arms

All fits use package defaults (75 trees, chi(1.5, 2), 500 burn-in; 4 chains, or 10 for the effect models,
as bartCause), `n.threads = 1`, and are read at several lengths from one run, a longer run's first draws
being a shorter one's (the binary interval study's prefix check).

- **shipped**: 207f740e with the move's toggle flipped off, the sampler bartcore ships. Readings: default
  (500 draws a chain after 500 sweeps of burn-in), equal time (more draws, the binary interval study's
  frozen choice, at the cell's time ratio), long (2,000 draws a chain, or 1,000 for 10-chain fits). Every
  replication.
- **move**: 207f740e as is, the same readings. In tranche 1 on reference replications only, where its
  default is the k-fast control from the same start; in tranche 2 on every replication of a cell that
  shows or is suspected of a problem (point estimate at or past the bar), where pricing needs its SBC
  readings.
- **oracle** (SBC cells): k fixed at the true k0, every replication. Given the data and k0, the true trees
  and leaves are exactly a draw from the posterior given k0, so this arm covers at nominal by construction
  once its trees mix: it is the "k known" sampler. D, the shipped arm's signed difference from it, is what
  is k's. A **long oracle**, 10 times every length, on 25 replications of S1 (factor covariates), S5 (250
  predictors), S6, S7 and E1 (missing values) in tranche 1, and of any other regime whose default oracle
  is half the bar or more from nominal: a long oracle off nominal is a mismatch between the prior draw and
  the fit (a harness defect, which stops that regime), a long oracle on nominal with a default one off is
  the trees' cost.
- **fixed k**: k fixed at the 10th, 50th and 90th percentiles of its posterior (the ruled arms, every
  chain), and at the 5th, 30th, 70th and 95th with one chain each, on every reference replication. With k
  fixed neither the probit step nor the move runs, so these fits are the same on both builds (a stage 1
  check) and run once. **Pooled stuck runs**: 200 random sets of as many chains as the user's fit has,
  drawn without replacement from the 16 (or 34) fixed-k chains with each chain's chance proportional to its
  quantile's share of k's posterior (the seven points stand for the bands 0 to 0.075, 0.075 to 0.2, 0.2 to
  0.4, 0.4 to 0.6, 0.6 to 0.8, 0.8 to 0.925, 0.925 to 1) divided by that quantile's chain count: a run whose
  chains each sit still at a k drawn from its posterior. Its mean mass is the stuck-run reading. The ruled
  three are also read alone. In a stationary 4-chain run some chain sits outside [q10, q90] with
  probability 0.59, which is why the bound rests on the pooled runs and not on the three.
- **stationary start**: the last 500 unthinned sweeps of each reference chain (4 chains, or 10, of each
  kernel), a default-length run of that kernel begun from its stationary distribution, at no extra fit.
- **D2 only**: a third pipeline with the move in the propensity stage and the shipped build in the outcome
  stage, to attribute a weighted-effect difference to the stage it comes from.
- **stan4bart arms** (E1, E2): stan4bart's own defaults (4 chains, 2,000 iterations, 1,000 warmup, one BART
  sweep per iteration) under the shipped and move libraries; k fixed through its BART arguments; equal
  time by total fit time, since it has no sampler-only timing.

**Where the k quantiles come from.** From the validated reference's pooled k draws, both kernels' chains
together, read only after the reference's final acceptance (after any escalation). That is not circular:
the arms under test are default-length runs, and a reference counts only when the shipped and move
kernels' long runs agree, k's quantiles included, so the quantiles do not rest on either kernel's being
right. When a reference fails, its quantiles come from its move half alone; those fixed-k readings are
flagged, enter only the 80 percent rule's sensitivity reading, and cannot make a cell "bounded". They never
come from the default fits under test. Under SBC the oracle needs no quantile.

**Equal running time.** The time ratio r (move over shipped, sampling time) is the binary interval study's
per cell for D1, D2 and B, from its timing pass (1.006 to 1.30). The SBC and embedded cells are timed the
same way in stage 1: a sampler-only run of 500 iterations, three alternating repetitions per build, on
three replications per regime, four jobs at once. The shipped equal-time reading has 1,000 r iterations a
chain.

## Truth construction

### SBC on real covariates

Per replication:
1. Take the regime's covariates (rotating over its datasets) and, where the regime's n is smaller, a fresh
   subsample; for effect regimes, the treatment and the shipped propensity (S7, S8), or a treatment drawn
   from a prior propensity (S9, below).
2. Build a dbarts sampler on that design with bart()'s defaults (75 trees, the default tree prior, k under
   chi(1.5, 2), offset 0, the default cut points), draw k0 from chi(1.5, 2), then trees and leaves from the
   prior at k0 ([`sbcReplication`](../../benchmarks/R/sbc.R)'s route), and evaluate f0 at every scored point
   (both treatment values for effect regimes).
3. Draw y ~ Bernoulli(F(f0)) with F the link.
4. Apply the regime's selection to (x, y). On rejection go back to step 2.
5. Fit every arm with bart() on (x, y) as a user would.

**Why selection keeps coverage nominal.** Under the joint distribution of steps 2 and 3, the truth given
the data is exactly a posterior draw: (k0, trees, leaves) given (x, y) follows p( . | x, y). Any acceptance
rule that reads only (x, y) leaves that conditional untouched, so given acceptance the truth is still a
posterior draw and an exact sampler's intervals cover at nominal, on average over accepted replications.
Selection changes which posteriors are probed, not their calibration. The same argument fixes three rules:
no acceptance test, row label or subgroup may read the truth; row labels (interior or tail) come from the
shipped arm's long reading, pooled with the move's where the move runs (never the oracle's, which knows
k0), a function of the data and the chains' own randomness, never from the truth; and coverage
conditional on the truth (for example, among truly null effects) is not scored under SBC. The oracle arm's
nominal coverage in every accepted set is the check.

**S9's two-stage truth.** The treatment z is drawn from a probit BART prior propensity e0 on the
covariates, with the S9 selection; then f0 from the outcome prior on (x, z) and y from f0. The joint
density factors as p(e0) p(z | e0) p(f0) p(y | f0, z), so the posterior given (z, y) is p(e0 | z) times
p(f0 | z, y): separate propensity and outcome fits give the exact joint posterior, and the weighted effects,
computed per draw from independent draws of each (and the data), are posterior quantities with nominal
coverage. A propensity entering the outcome model as a covariate would break this, which is why S9 leaves
it out and D2 carries it.

**stan4bart's prior** (E1, E2): BART from the dbarts prior as above, group intercepts from a normal with
the group standard deviation drawn from stan4bart's prior on it, and the parametric intercept from its
prior, each transcribed from stan4bart's prior code. Stage 1 checks the transcription against a stan4bart
fit with every weight zero, which samples the prior.

**E3's truth** is the S6 construction; 30 percent of rows are then masked at random and imputed by the
loop.

### Near 0 and 1

The binary interval study's truth near 0 and 1 was a posterior mean, which on near-separated data lies
outside its own posterior's interval, so any sampler that reached the posterior missed it. Here:
- the SBC truth is a draw from the posterior, typical of it by construction, so that cannot happen;
- intervals and coverage are computed on the latent scale and every stored draw is latent, so no
  probability rounds to 0 or 1;
- the interval score on the probability scale, which weighs a miss by its size, is reported beside
  coverage;
- for designed truths (D1, D2) the truth is a designed surface, not a posterior mean, and only secondary
  readings use it; the widened coverage (0.001 on the probability scale) is reported beside the plain one
  for rows near 0 or 1.

### References

- **Run**: 8 chains of each kernel (10 for 10-chain cells), each 5,000 sweeps of burn-in, 20,000 sweeps
  recorded at thin 20 (1,000 saved draws a chain), then 500 sweeps saved unthinned (the stationary tails),
  through the sampler API so a chain can be continued. Latent-scale single-precision draws of the scored
  quantities, k (and the group standard deviation) are written. stan4bart has no continuation: its
  references are separate long stan4bart runs (5,000 warmup, 20,500 sampling iterations), every iteration's
  scored quantities recorded through its callback, thinned afterwards, the last 500 iterations being the
  tail; an escalation is a longer rerun of the shipped half.
- **Acceptance** (the orchestrator's reading of dec-B462's "agree and converge", below): the two kernels'
  halves agree on every scored quantity's mean and 5th and 95th percentiles, and on k's 5th, 10th, 30th,
  50th, 70th, 90th and 95th percentiles, within 4.5 combined Monte Carlo standard errors; rank-normalized
  split-Rhat over all chains at most 1.01 for every scalar estimand and at most 1.02 at the 99th percentile
  of row quantities; tail effective draws at least 400 for every scalar estimand and at the 1st percentile
  of row quantities (the thresholds of Vehtari et al. 2021; the 1st percentile replaces the pilot's minimum
  over rows, which its review judged strict); k's effective draws at least 100, the precision that places
  its 5th and 95th percentiles within about 2 percentile points, which the fixed-k arms need. k's own
  split-Rhat is reported and not required.
- **Escalation, only where it can pass** (the trigger amended in [Pre-registration](#pre-registration)):
  where the halves disagree on k or k's effective draws fall short, the shipped chains continue to 100,000
  recorded sweeps when the first pass's effective draws per sweep, for k and for the failing quantities,
  project a pass at that length; otherwise the reference fails at once. The shipped kernel's k needs 5,000
  to 30,000 sweeps per independent draw near separation
  ([k-mixing-pg-families.md](../design/k-mixing-pg-families.md)), so 8 chains of 100,000 give about 27 to
  160 effective draws of k, which can meet 100 and could never have met the pilot's 400. A reference
  failing after escalation does not count, and the failure is recorded as a finding (the shipped kernel
  does not reach its own long run there).
- **Storage and order** (the audit set amended in [Pre-registration](#pre-registration)): full draws are
  kept until every reading that depends on the reference (band mass and other per-draw functionals need
  joint draws) has been scored, and every fixed-k job waits for the final acceptance. Then kept: per
  quantity, 1,000 quantiles of the pooled draws (single precision), k and scalar estimands in full, and
  full draws for 20 references chosen in advance (audit); about 1.2 MB a reference.
- **Check on the references as a class**: under SBC, the references' own coverage of the truth must be
  nominal within noise over a regime's reference replications; a regime where it is not shows acceptance
  passing references it should not, and its reference readings are withdrawn.

## Replications and power

Per-replication spreads below are planning figures; tranche 1 measures them and tranche 2's counts follow
from them ([Order of work](#order-of-work)).

| reading | assumed spread a replication | tranche 1 | full |
|---|---|---|---|
| SBC, row or unit classes (f, individual effects): D | 0.08 (the pilot's paired spread was 0.05 between two drawn-k arms; against the oracle it is larger) | SE 0.008 at 100: problem from about 0.029; D inside the bar when the estimate is within 0.004 of 0 | SE 0.004 at 400: from 0.02; within 0.012 of 0 |
| SBC, scalar classes (average effects, subgroups, ratios, group sd, weighted effects): D | paired discordance 0.10 | SE 0.032: problem from about 0.11 | SE 0.016 to 0.018: from about 0.06; no problem not resolvable by coverage alone |
| reference readings: S, stuck runs, transients | 0.03 where nothing is wrong (the corrected pilot reading sits at nominal with little spread); 0.15 where something is | 12 references: SE 0.009 at t 2.20, so "no problem" needs an estimate near 0; problem from about 0.04 at t of about 5 | 24: SE 0.006, no problem below 0.0075 (power about 0.90 at a true 0); problem from about 0.027 where the spread is 0.03, 0.13 where it is 0.15 |
| D1 contrasts, by reference | as above, 16 references a link | | as the 24-reference column, a little wider |

So scalar estimands are resolved for "no problem" by references, not by SBC coverage, which for them can
show only large problems; and tranche 1's 12 references can show a reference problem but almost never
clear one, which is why every unresolved tier 1 and 2 cell expands. That is a property of one yes-or-no per
replication and of small reference counts, not a choice of bar.

## Cost

Per-fit seconds come from the binary interval study's pilot (its per-fit times, plus its scoring, half of a
small pair job); reference seconds are sweeps times the timing pass's per-sweep times, which the pilot's 52
references matched. Escalation is priced on a quarter of the references in the near-separated and
small-sample regimes and a sixth elsewhere, as the predicted-pass rule should leave. Core-hours on the box:

| part | tranche 1 | tranche 2, expected | tranche 2, at most |
|---|---|---|---|
| stage 0, the screen on stored outputs | 0.2 | | |
| stage 1: builds, checks, SBC self-checks, timing, smoke | 4.0 | | |
| long oracle (S1, S5, S6, S7, E1; others if needed) | 3.0 | 2.3 | 4.7 |
| rare-outcome co-selection (draws and tree scores) | 1.0 | | |
| S1 near separation (tier 3: does not expand) | 2.1 | 0 | 5.4 |
| S2 rare, 500 | 2.0 | 4.1 | 5.3 |
| S3 rare, 10,000 | 7.1 | 4.2 | 5.5 |
| S4 small n, strong | 1.1 | 2.4 | 3.1 |
| S5 many predictors | 6.0 | 5.4 | 6.9 |
| S6 control | 1.8 | 4.1 | 5.3 |
| S7 effects | 5.6 | 8.2 | 10.0 |
| S8 rare effects | 8.9 | 8.6 | 11.4 |
| S9 propensity, weights (probit) | 3.5 | 5.1 | 6.2 |
| E1, E2 stan4bart (probit) | 6.6 | 5.6 | 6.8 |
| E3 dbartsSampler loop | 2.3 | 2.2 | 2.9 |
| D1 multiple treatments, by reference | 20.4 | | |
| D2 weighted effect, by reference (probit) | 4.4 | | |
| B real outcomes | 8.6 | | |
| total | 88.4 | 52.2 | 73.5 |

**Expected total about 141 core-hours**, about 40 hours of the box at four jobs (15 percent added for
shared memory bandwidth). Expected tranche 2 assumes every tier 1 and 2 cell expands, under both links
(dec-B463), as the power table says they will, with the move run on every replication in a third of them,
and the long oracle extended to half the remaining regimes. At most about 162, or 167 with the colleagues'
case. **The stop is at 211**, 1.5 times the expected total, and covers every expansion. Under the stop,
cuts go in this order: S3's references, tranche 2's references in tier 2 cells, tranche 2's SBC
replications in tier 2 cells.

Against the critique's estimate of about 165 after its changes: the difference is the move run only on
reference replications until a cell shows a problem (the critique's optional cut, about 8 in each
tranche), D2 as a probit reference cell (about 5) and S9 probit only (about 11 over both tranches).

Of tranche 1, references are about 55 core-hours. Disk: about 3 GB of new output (fit summaries about
1 GB, compacted references about 0.5 GB, the 20 audit references about 1 GB), with full reference draws of
up to about 2.5 GB held at a time until they are scored. No default-reading draws are stored: every
reading keeps its interval ends at four levels, mean, standard deviation and latent median per quantity,
which is all the pointwise mass reading needs; bands are scored before the reference is compacted. The box
has about 21 GB free with the binary interval study's 18 GB kept until its write-up. Memory: S3's fits
hold 4 x 2,000 draws of 10,000 fitted values (0.64 GB) and D1's as the pilot's; both run at most two at a
time.

### Considered: a lean variant

A lean variant of about 88 core-hours was considered and not taken (the orchestrator's call, under the
maintainer's standing permission to raise compute when it gets a better answer): logistic dropped from
S7, S8, E3 and D1, and tranche 2 only for tier 1 cells. It would leave tier 2 cells at tranche 1
resolution, so their scalar estimands (ratio estimands at rare rates above all) could show only a large
problem and never be cleared, and it would say nothing on logistic effects, logistic embedded use or the
logistic multiple-treatments design, which is where k mixes worst (dec-B463). The full design resolves
every tier 1 and 2 cell, under both links, to problem or no problem at the bar wherever references pass;
near separation, tier 3, is read at tranche 1 resolution (a problem of about 0.03 in row classes and 0.04
in reference readings can be shown; otherwise it mostly ends unresolved).

## Order of work

0. **The screen, stage 0** (no fits; minutes on one core of the box). On the binary interval study's 52
   pilot references, which store k with every draw, and its corrected re-read: for every scored quantity,
   the share of its posterior variance explained by log k (draws binned by k decile), and its conditional
   mean and 90 percent interval at the bins holding k's 10th and 90th percentiles, in posterior standard
   deviations and as reference mass. Interior rows and effects only, as the re-read defines them. Output: a
   per-dataset table of the bound's proxy and a stated prediction, per regime and estimand class, of which
   cells tranche 1 will find bounded. It gates nothing; it is written down before tranche 1 is read so the
   main run can be checked against it. Caveat stated with it: within a reference, conditioning on k only
   approximates a fixed-k fit, and 34 of the 52 references failed the pilot's acceptance.
1. **Stage 1** on the box ([Gates before the main run](#gates-before-the-main-run)): builds, checks,
   timing, the selection bands and the cost per accepted replication, the SBC null and discrimination runs,
   the stan4bart prior check, the smoke.
2. **Review** of the harness and smoke; fixes.
3. **Tranche 1**, S6 and S1 first (the two ends), then the rest in random cell order. Job dependencies: a
   reference's escalation decision, then its final acceptance, then its fixed-k and move fits, then their
   scoring, then compaction.
4. **The tranche 1 reading**, at the levels of [Bounds](#bounds), written up and reviewed. A cell expands to
   its full counts when it is tier 1 or 2 and unresolved, or when it shows or is suspected of a problem
   (pricing needs its size and transient share to two standard errors of 0.01); tier 3 cells stay as read.
   Full counts are capped at the table's and lowered where the measured spread resolves the bar at fewer
   ((z_a + z_b)^2 s^2 / 0.02^2 with z_a 1.96 and z_b 0.84). The move joins every replication of a cell
   that shows or is suspected of a problem.
5. **Tranche 2**, tier 1 cells first.
6. **The final reading**, the map, the where-it-matters analysis below, the finding, a ledger entry for the
   calls made running it, and the result to the maintainer.

**Where it matters.** Beyond the cells, across all SBC replications: the per-replication D (and, on
references, S) against data-only features (n, outcome rate, events, separation score, predictors, the
shipped run's own k split-Rhat) in an additive model, exploratory. It says whether a data-defined trigger
could find the fits that need fast k moves, which is a pricing option (below).

## Gates before the main run

Harness checks, each a stop if it fails:
- **Builds**: both from 207f740e as the binary interval study made them (its box libraries reused if they
  pass, else rebuilt with `--preclean`); stan4bart's bartcore branch installed against each.
- **The shipped arm is what ships**: the shipped build bit for bit equal to the bartcore tip of the day on
  the study's paths (both links, numeric and factor predictors, missing values, test rows, 4 and 10
  chains, the sampler API's one-sweep loop with setResponse, a stan4bart probit fit). The study runs after
  the merge to main, so the engine may have moved since 207f740e's base; if draws differ, the orchestrator
  rebases the move's branch onto the tip and both arms are rebuilt from it.
- **The move differs**, and **fixed-k fits are equal on both builds**, bit for bit; the prefix property on
  both.
- **SBC self-consistency** on every regime's design: the SBC harness's fit-consistency and k hand-back
  checks extended to logistic and to the real designs (factors, missing values, 250 predictors, effect
  designs at both treatment values); the prior-draw sampler's cut points equal to bart()'s on the same
  data.
- **SBC null run**: S6 at 100 replications, both links, the oracle and the long oracle within noise of
  nominal coverage. **Discrimination runs**, each of which must flag: the outcome drawn through the other
  link, read on the oracle (whose fixed k cannot absorb the link's scale); and the truth drawn from a prior
  of 50 trees and fit at 75, read on the drawn-k arms, a mismatch in the function's shape that a drawn k
  cannot absorb. A null failing is a harness defect, not a finding.
- **stan4bart prior check**: 2,000 draws of a zero-weight stan4bart fit against the harness's prior draws
  of the group standard deviation, intercept and BART fit at five quantiles each.
- **Selection**: the bands set as [Regimes](#regimes) says, and each regime's cost per accepted
  replication measured against its fit cost, before any arm runs.
- **The pilot's stored readings**: that the binary interval study's pilot stored D1's shipped default
  reading for replications 1 to 10 (if only its equal-time reading was kept, the secondary reading uses
  that and says so).
- **The binary interval study's defects, not repeated** (its review, items M1, S1 to S4 and the NITs):
  every job list covers every replication of every phase, propensities included, and a missing input
  stops the job rather than dropping a pair; no probability is stored or scored in single precision;
  every summary is computed per column with its own missing-value handling, and a check asserts that each
  summary's row count is the planned count; multiple-treatments references are scored on the six
  contrasts from the four arms' latent draws; no reading writes to another reading's path; jobs are
  written in tranche order so no tranche runs past its reading, and no job reads a reference before its
  final acceptance; a zero spread is not treated as no noise; the bootstrap resamples whole replications
  across cells.
- **Truth-blind labels**: a check that no acceptance rule, row label or subgroup reads f0, k0 or e0.
- **Scoring halves**: a check that every kernel's readings are scored against the other kernel's half
  and a tail never against its own chain.

Then a **smoke**: every cell, one replication, both links and builds, at a tenth of every length, on the
box, through every reading. Then **review**, and only then tranche 1. A pre-registration commit of this
plan's rule, bands and bars precedes the first study fit.

## What the result leads to

**A problem shown** (k's, in some cell) opens the pricing step, a separate plan, which reads this study's
outputs and measures only variants:
- per problem cell, the effective draws per second of the affected estimands under both kernels (from the
  long readings and the stationary tails), and the breakdown of the shortfall;
- candidates: the built move every fourth sweep (207f740e); every sweep; during burn-in only, which removes
  the transient and nothing else, so it suffices exactly where k's share of the transient is the whole
  problem; a longer default burn-in; the move gated by a data-defined trigger (outcome rate or separation
  score), exact because it reads the data alone, as the collapsed plan's size switch does; the collapsed
  plan's size and frequency options; an opt-in argument;
- what follows by tier: a tier 1 problem prices a default (the cheapest candidate that clears the bar
  where the problem is, at its cost on every fit); tier 2 alone prices a scoped default (trigger) against
  an opt-in, by the share of fits in those regimes, which the maintainer estimates; tier 3 alone leads to
  documentation and a saturation-keyed diagnostic, with the move at most an opt-in;
- the costs beside it: the collapsed plan's 1.4k to 2.9k lines, a posterior-changing slice, consumers'
  expectations re-recorded, and, for the probit half, every binary stan4bart fit.

**No problem** (every tier 1 and 2 cell no problem, bounded or measured):
- written down in docs/design/binary-k-mixing-need.md: dec-B420 settles that k is drawn, not that it must
  mix fast; the map; the stuck-run bound per estimand class, which holds for any sampler; that the reported
  k and probabilities near 0 or 1 depend on chain length near separation; the oracle and tail readings as
  the trees' own cost;
- the help pages for bart() and the binary families say the reported k and extreme probabilities near
  separation depend on chain length;
- the move stays parked on wt/collapsed-move-plan at 207f740e, its revival cost (a rebase and its gates)
  stated; the maintainer is asked whether dec-B404's and dec-B420's "built before 1.0-0" stands;
- reopening triggers, in TODO: a reproducible user fit whose reported interval or effect moves with chain
  length or seed by more than the bar; a change to k's default prior or a new binary family; an embedded
  use with fewer BART sweeps per outer iteration than stan4bart's one, or a new consumer embedding BART; the
  tree-mixing proposals landing (they change the comparator); a request to report k as a parameter; a fit
  outside the tested envelope (rare outcomes beyond 10,000 rows, more than 250 predictors, fewer than 8
  groups).

**No problem shown, some tier 1 cells unresolved** (pre-stated, as the likeliest outcome for some scalar
estimands): the default action is the no-problem path above. The maintainer asked first to confirm that a
problem exists; an unresolved cell has not confirmed one. Two additions: each unresolved tier 1 cell goes
to the maintainer with its estimates, bounds, the reason it is unresolved and the cost to resolve it; and a
cell whose point estimate is at or past the bar (suspected) gets a targeted extension within the stop
before the write-up, after which it is read again by the same rule. Unresolved cells are added to the
reopening triggers by name.

**Unresolved cells in tiers 2 and 3** are listed in the finding with the counts each would need.

The study's diagnosis (gaps the oracle shares, references that fail even after escalation, the trees'
cost) goes to TODO tree-mixing-proposals whatever the reading.

## Harness

benchmarks/R/binary-k-mixing-need.R, tracked: the regime registry and seeds, the SBC construction on real
designs and its selection, the arms and their readings, the pooled stuck runs, the reference runner with
predicted escalation and tails (and stan4bart's callback recording), scoring by halves, the readings of
[What counts as a problem](#what-counts-as-a-problem), the tranche rule, the stan4bart and dbartsSampler
drivers, and actions `screen`, `check`, `timing`, `jobs`, `fit`, `reference`, `read`. It reuses the binary
interval study's generators (A6, ACIC 2016 settings, IHDP, Twins), loaders and measures, moved into a file
both harnesses source once that study's harness is on bartcore, as binary-datasets.R was; the SBC harness's
prior draw and checks; and [benchmarks/R/uci-binary.R](../../benchmarks/R/uci-binary.R). The box driver, in
scratch, is the binary interval study's queue with job dependencies. Outputs stay in scratch.

## Calls made in planning

**The orchestrator's**, on the critique's contested points, 2026-10-17:
- **dec-B462's "converge" is read on the scored quantities, plus k's quantiles at the precision the
  fixed-k arms need (about 100 effective draws), not on k's own split-Rhat at 1.01.** Requiring k itself
  to converge selects the replications where the shipped kernel already mixes, which biases the reading
  toward no problem. The ruling's other guard stays: both kernels' long runs must agree, k's quantiles
  included.
- **A cell's reference readings count only when most of its attempted references pass**: 80 percent, the
  critique's figure, with the added condition that the reading survive scoring each failing replication
  against its move half. Otherwise the cell is unresolved, with the reason stated. A fixed share alone
  cannot rule out that the failures are the hard replications; the sensitivity reading tests that, and 80
  percent keeps it from carrying more than a fifth of the cell.
- **The designed cells restructured as the critique proposed**: D3 cut; D1 read by references with the
  fixed-k and transient readings, its pilot replications a secondary reading; a slot for the colleagues'
  exact case (D1x). D2 is made a reference cell by the same reasoning.
- **A default action pre-stated** for "no problem shown, some tier 1 cells unresolved"
  ([What the result leads to](#what-the-result-leads-to)).
- **S9's logistic half dropped**, and D2's with it (the propensity stage bartCause runs is probit).
- **The full design, not the lean one**, under the maintainer's standing permission to raise compute when
  it gets a better answer ([Considered: a lean variant](#considered-a-lean-variant)).

**The maintainer's**: in the use tiers a logistic cell counts as the same cell under probit (dec-B463);
the plan had proposed counting it one tier rarer.

**The critique's other items, taken**: the nominal comparator and the breakdown of the shortfall (its M1);
escalation only where it can pass (M2); the rare regimes co-selected for weak signal, reported by
separation score and gated on cost per accepted replication (M4); the pooled stuck runs over seven
quantiles beside the ruled three, and the bounded reading's claim narrowed (S1); G and signed D, pooling
over replications, the bootstrap variance in the noise correction (S2); k's share of the transient against
the move from the same start (S3); the long oracle (S4); D1's bands in tier 1 (S6); the inverse-weighted
estimates in S9 (S7); bounds from the bootstrap standard error, t quantiles for reference counts (S8);
stan4bart references through its callback and separate long runs (S9); full reference draws kept until
every dependent reading is scored, fixed-k jobs after final acceptance (S10); the cost re-estimated (S11);
the default action (S12); and the NITs: sweeps stated, the tests counted (z = 3.60), tails and defaults
paired on one half, the pilot's single precision acknowledged and its stored reading checked, Twins' rows
as stan4bart keeps them, a tree-count mismatch added to the discrimination runs, B tiered per dataset, the
move's k split-Rhat quoted with the panel's 1.6, S5's dense truth stated.

**Taken in part, with the critique's alternative**:
- The critique judged the link-swap discrimination run absorbed by a drawn k and proposed a tree-count
  mismatch instead. Kept both: the link swap is read on the oracle, whose fixed k cannot absorb it, and the
  tree-count mismatch on the drawn-k arms.
- The critique offered Holm's procedure as an alternative to counting the tests. Kept Bonferroni at z =
  3.60: simpler to state in advance; Holm's procedure starts at the same level and gains only after a
  first problem is shown, when the reading is no longer in doubt.
- The critique proposed q05, q50, q95 in place of the ruled quantiles. Kept the ruled q10, q50, q90 as
  arms and added q05, q30, q70 and q95 at one chain each for the pooled runs, which carry the bound.

**Mine**, each with its alternative, open to review:
- **Problem readings on the shipped default**, not at equal time: users run defaults; equal time is a
  pricing comparison. Alternative: read at equal time, as the binary interval study did.
- **The oracle (k at its truth) as the SBC comparator** beside nominal. Alternative: nominal only, which
  folds the trees' cost into k's.
- **Stationary start from reference tails** scored against the other kernel's half. Alternative: dedicated
  fits from saved reference states, about 6 more core-hours, independent of the reference's chains.
- **k quantiles from the validated reference**, the move half flagged as fallback. Alternatives: a k grid of
  fixed-k fits with the marginal likelihood ratio across neighbors estimated from each draw's leaf count and
  leaf sum of squares, which uses no drawn-k chain at all but is a new estimator needing its own validation;
  or no fixed-k reading where no reference passes.
- **The move on reference replications only in tranche 1** (the critique's optional cut). Alternative: on
  every replication from the start, about 8 core-hours more a tranche, giving the move's SBC readings in
  every cell.
- **Two tranches with a pre-registered expansion rule**. Alternative: full counts everywhere, about 162
  core-hours, with no look.
- **One bar, 0.02, for every class, bands included.** Alternative: a wider bar for bands and scalar
  estimands where the counts cannot resolve 0.02.
- **Default offset 0 with rare outcomes by selection.** Alternative: a fixed nonzero offset in prior and
  fit, which is not the user's default model but would reach weak-signal rare outcomes without rejection.
- **Large n at 10,000 rows.** Alternative: 20,000 (S3 about twice its cost) or 50,000.
- **E3, the dbartsSampler loop**, added beyond the brief because swapped responses are dbarts's other
  distinguishing use and data augmentation slows mixing; about 5 core-hours. Alternative: drop it.
- **No default-reading draws stored.** Alternative: store them for later re-scoring, about 30 GB, which does
  not fit the box.
- **Reuse of the binary interval study's pilot**: D1's replications 1 to 10, its references as the screen's
  input, its time ratios and its frozen equal-time choice.

## Pre-registration

Fixed 2026-10-29, after stage 0, stage 1 and the smoke and before any study fit. The harness is
benchmarks/R/binary-k-mixing-need.R as of the commit that adds this section; the bands below are written
into it.

**Unchanged by stage 1**: the bar, 0.02 for every class
([Bars and where they come from](#bars-and-where-they-come-from)); the
bounds, z 3.60 for a problem and 1.96 for no problem, Student's t at the same levels with one degree of
freedom fewer than a cell's passing references ([Bounds](#bounds)); the reading of a cell, with its 80
percent pass rule and the reading against each failing replication's move half
([The reading of a cell](#the-reading-of-a-cell)); the tiers
([Reading by regime and weighing by use](#reading-by-regime-and-weighing-by-use)); reference acceptance and
escalation ([References](#references)), the escalation's trigger amended below; the expansion rule
([Order of work](#order-of-work)); the stop at 211 core-hours and the order of cuts under it
([Cost](#cost)).

**Selection bands**, the separation score's central 80 percent over the k study's simulated cells at an
outcome rate of 0.5 whose truth-based sign agreement is within 0.03 of 0.76 (9 cells) and of 0.90 (5
cells):

| band | separation score |
|---|---|
| moderate | 0.574 to 0.807 |
| strong | above 0.807, below 0.853 |
| near-separated | 0.853 to 0.915 |
| deterministic, rejected | 0.995 or more |

A real dataset of B scoring above 0.915 takes the near-separated tier.

**Cost per accepted replication**, on three replications a cell and link: every regime's selection costs
less than one of its fits, so every regime runs at its planned counts.

| regime | draws per accepted design, probit / logistic | seconds to select one | seconds of one shipped fit |
|---|---|---|---|
| S1 | 8 / 13 | 0.4 / 0.5 | 1.4 / 2.6 |
| S2 | 13 / 9 | 0.2 / 0.3 | 1.7 / 2.9 |
| S3 | 11 / 7 | 0.6 / 0.5 | 15 / 38 |
| S4 | 25 / 12 | 0.4 / 0.3 | 0.7 / 1.1 |
| S5 | 1 / 1 | 0.2 / 0.2 | 2.4 / 4.5 |
| S6 | 2 / 2 | 0.2 / 0.2 | 1.1 / 2.2 |
| S7 | 1 / 1 | 1.7 / 0.7 | 11 / 14 |
| S8 | 6 / 13 | 2.4 / 1.2 | 13 / 16 |
| S9 | 1 | 0.2 | 13 |
| E1 | 1 | 0.8 | 5.9 |
| E2 | 25 | 0.9 | 5.7 |
| E3 | 3 / 1 | 0.2 / 0.2 | 4.4 / 6.1 |

The near-separated band fills without deterministic draws (S1, 8 and 13 draws an acceptance).

**Time ratios**, move over shipped, the median of nine paired timings a cell and link: SBC and embedded
cells 0.98 to 1.29, mostly 1.04 to 1.14; D1, D2 and B take the binary interval study's (1.006 to 1.24).
The plan assumed 1.006 to 1.30.

**Projected cost**: tranche 1 about 80 core-hours from the timing pass's per-sweep times (58 of sampling,
9 of them escalation, plus 15 percent for shared memory bandwidth and 20 percent for scoring and
start-up), and about 77 from the smoke's job times scaled to full length (each job ten times its smoke
time less start-up, acceptance 3.8 times as its diagnostics scale, plus the escalation's 9). The plan has
84.2 for the same parts. Stage 0, stage 1 and the smoke used about 5.5 core-hours against the plan's 4.2,
the excess the second discrimination run and the smoke's acceptance passes. At the plan's ratio of the
whole study to tranche 1 (141 to 88), the expected total is about 130 core-hours, under the stop at 211.
Amended below.

**Gates**, all passed ([Gates before the main run](#gates-before-the-main-run)):
- Builds: the shipped build bit for bit equal to main on 22 paths (both links; numeric, factor and missing
  predictors; test rows; 4 and 10 chains; offset; fixed k; the sampler loop; stan4bart probit and fixed
  k); the move differs on every drawn-k path; fixed-k fits equal on both builds; the prefix property on
  both, stan4bart included.
- SBC self-consistency on all 21 regime and link designs: cut points, fits at train and test rows, k
  handed back, k moving.
- Null run, S6 at 100 replications a link: the oracle and long oracle cover interior f at 0.900 and 0.901
  (probit), 0.891 and 0.893 (logistic). The logistic default oracle sits 2.5 standard errors below
  nominal, inside the call below; a logistic G carries that offset, which D, paired with the oracle, does
  not.
- Discrimination: as the call below, amended below.
- stan4bart prior check: twenty quantiles (the BART fit, the group standard deviation, two groups'
  intercepts) within 1.6 standard errors of the harness's prior draws.
- Selection: the bands and costs above. The pilot's stored readings: D1's shipped default reading for
  replications 1 to 10, under both links. The binary interval study's defects, truth-blind labels and
  scoring halves: checked, the last on the smoke.
- The smoke: all 1,655 jobs, every cell with one replication under both links and builds at a tenth of
  every length, ran through every reading; the compaction added after it then ran on its 81 references
  (5.1 GB to 0.9 GB, the audit references kept). It escalated no reference, so escalation was forced on
  an S6 and an S9 reference and now runs through acceptance, fixed k, scoring and compaction. Its single
  replication read eleven SBC cells as problems on the binomial floor of the standard error alone, which
  led to the minimum below.

**Calls made implementing it**, each with its alternative:
- Bands from the k study's cells at an outcome rate of 0.5 only. Alternative: every rate, whose scores
  overlap and give no clean band.
- Bands fixed at the bench box's values; the same computation elsewhere can differ in the third decimal.
  Alternative: recompute where the study runs.
- E3 draws its missing-outcome mask before the truth, and selection reads the observed rows only.
  Alternative: all rows, which reads outcomes the loop imputes.
- D2: each pipeline fits its own propensity, the reference uses the shipped pipeline's, and the mixed
  pipeline is the move's propensity with the shipped outcome model. Alternative: a reference per pipeline,
  twice the cost. Amended below.
- D2's TMLE estimate is not read: it is not a posterior quantity. Alternative: its interval's reference
  mass.
- B: held-out f at the first 100 test rows of the k study's split, and its six seeds as its replications
  against its one reference. Alternative: all test rows, as the pilot; pooling problems by regime.
- E2's states pool into the data's nine census regions. Alternative: merge two of them to make eight.
- stan4bart references keep their fits instead of recording through the callback, since a callback fit
  keeps no k. No alternative records k.
- The stan4bart prior check fits one row at weight 1, all others at 0, and weights the harness's prior
  draws by that row's likelihood: stan4bart's probit takes weights of 0 and 1 only. Alternative: a tiny
  weight, which it refuses.
- SBC row labels average the shipped and move long readings' 5 and 95 percent points where the move runs.
  Alternative: pooled draws.
- Under SBC, band classes are read on D alone: a band built from a draw count is not nominal by
  construction. Alternative: G as for other classes.
- Mass readings sum the two halves' counts and never build a pooled copy. Alternative: a pooled copy,
  about 5 GB of memory for a D1 reference. Amended below.
- In 10-chain cells, a tail's reference without its own chain is the other kernel's half alone.
  Alternative: the pool less one chain.
- Escalation is not built for stan4bart and dbartsSampler-loop references: one that fails acceptance stays
  failed and counts against the 80 percent rule. Alternative: a longer rerun of its halves. Amended below.
- The null run passes when every reading is within 2.81 standard errors of nominal, a 0.005 two-sided level
  over its twelve readings. Alternative: 1.96 a reading. Amended below.
- A discrimination run flags when the 95 percent bound of interior f's coverage excludes 0.90. The plan's
  two runs did not both flag under both links: the outcome through the other link flags under probit
  (coverage 0.70) but not logistic (0.896), whose oracle at the true k shrinks its fit back onto the
  truth's scale; the 50-tree truth flags under neither, its coverage equal to the null run's (0.899 and
  0.884 against 0.902 and 0.885), since the leaf prior sets f's prior variance by k whatever the number
  of trees, and a 75-tree fit absorbs the shape. The gate is instead a truth whose f is doubled while the
  oracle keeps the drawn k, a wrong k by a factor of 2, read on the oracle under both links (coverage 0.72
  and 0.71, both flagged); the plan's two runs are read and reported. On the shipped arm the drawn k
  absorbs most of the doubling (0.88 and 0.86, against the null run's 0.90 and 0.88). Alternative: stop
  for a ruling on the gate. Amended below.
- Tranche 1 runs in a fresh run directory: the null run's S6 fits are gate fits, not study fits.
  Alternative: keep them as S6's first replications.
- stan4bart and the loop are timed by total fit time, having no sampler-only path.
- The smoke reads stage 1's measured time ratios. Alternative: a fixed ratio, which skips every cell's
  lookup.
- An SBC cell reads as unresolved on fewer than 20 replications: a bootstrap over so few cannot see the
  spread between replications. Alternative: the binomial floor alone, which read the smoke's single
  replication as problems.
- A reference is compacted once scored: per row quantity 1,000 quantiles of the pooled draws, scalar
  estimands and k in full, band classes (joint draws of the row quantities, already scored) dropped, saved
  samplers removed; 20 references chosen in advance keep their draws. Alternative: keep the samplers of
  the audit references, so their chains could be continued, about 40 MB a part. Amended below.
- After an escalation the shipped chains are five times the move's length; the joint split-Rhat and tail
  effective draws take them thinned to the move's length over their whole run, while the halves'
  agreement uses every draw. Alternative: the halves' effective draws summed, with no joint Rhat.
  Amended below.
- A reference that does not escalate takes its first acceptance as final. Alternative: recompute it, about
  9 core-hours of tranche 1.
- D1's secondary reading on the pilot's replications 1 to 10 is checked present but not yet built into the
  reading; it is built before tranche 1 is read. Amended below.

**Amendment, 2026-10-29, before tranche 1**, after the review of the pre-registration and before any
study fit: each change, what it replaces, and why. The harness is benchmarks/R/binary-k-mixing-need.R as of
the commit that adds this amendment.

- **Discrimination gate.** The link run is the gate again. The doubled-truth gate that replaced it rested
  on a scoring bug: the link run scored the oracle's interval, on the fit link's latent scale, against the
  truth's latent value on the other link's scale, where it should have scored it against the truth's
  probability carried onto the fit link's scale. Rescored, the oracle's interior f covers at 0.944 under
  probit (se 0.003) and 0.761 under logistic (se 0.006): both flag. (The review's 0.946 and 0.734 label
  interior rows by the oracle's own interval; the harness labels them by the shipped arm's.) The doubled
  truth read on the oracle (0.724 and 0.706, se 0.007 and 0.006) stays as a plumbing check that must flag
  too, its shipped arm (0.879 and 0.865) reported as a contrast. The 50-tree run is dropped: the leaf prior
  sets f's prior variance through k whatever the number of trees, so it cannot flag.
- **Graded controls**, a new gate: on the null run's S6 designs, the oracle at k fixed at 1.25 and at 0.8
  times the true k, read as the shipped arm is read, must read "problem" in interior f under both links. A
  problem reading needs both G and D at or past the bar of 0.02 and each more than 3.6 standard errors
  from zero. It passes, the margins over the bar in standard errors: at 1.25 (probit, logistic) G is 0.042
  and 0.057 (se 0.004 and 0.005; margins 5.2 and 7.6) and D -0.042 and -0.047 (se 0.002; margins 13 and
  15); at 0.8, G is 0.024 and 0.021 (se 0.003; margins 1.4 and 0.4) and D +0.025 and +0.031 (se 0.0014 and
  0.0016; margins 3.3 and 6.8). The 0.8 control passes thinly under logistic: its G clears the bar by 0.4
  standard errors, which follows from the rule rather than from a defect. Measured by the review, D clears
  the bar reliably from about 1.5 times the bar, while a reading of over-coverage near the bar is close to
  a coin flip; scalar reference classes have no control and resolve only from about 0.045 at 12
  references.
- **Reference-path control**, a new gate: on tranche 1's count of S6 and S7 probit references (12 each,
  full length, fresh seeds), every chain of the user's fit at k fixed at 1.25 times the reference's median
  k, read through the reference reading as the shipped default is, must read "problem" in S6's interior f
  and S7's unit effects. It passes: the control's S is 0.052 (se 0.0015) in S6's interior f and 0.032 (se
  0.004) in S7's unit effects, against the shipped default's 0.002 and 0.002, both bounded. 11 of 12
  references passed acceptance in each cell; the two that did not escalated and failed again, one on a
  single row's halves (z 4.9 of 600 quantities, its effective draws over 15,000), one on rows' Rhat (1.04
  at the 99th percentile). It used 1.6 core-hours. The control also runs on tranche 1's S6 and S7
  references and is reported.
- **Null bound.** The null run's bound is Bonferroni over its 12 readings at a two-sided 0.05: |z| at most
  2.87, in place of 2.81. The run passes (largest |z| 2.52, the logistic default oracle). Both bounds,
  2.81 and 2.87, were set after the logistic -2.52 was seen, so neither settles it: the review's null run
  of 600 further replications does, putting both links' oracle at 0.898 (se 0.0015), and its longer oracle
  shows the default oracle's finite length costs about 0.2 points of coverage: the logistic offset is
  noise and length, not a defect.
- **Seeds.** Every seed drawn while an SBC truth is made, and every arm's seed, carries the run's tag (its
  directory name; stage 1 has none), so tranche 1's S6 replications are fresh against the null run's gate
  designs. Designed and real designs keep the binary interval study's seeds.
- **Escalation.** A failing first pass escalates when the shipped chains' effective draws project a pass
  at the escalated length, whatever else failed; the halves' agreement, k's agreement and Rhat are judged
  at the final pass. [References](#references) escalates a disagreement on k or a shortfall of k's
  effective draws where the effective draws project a pass, but the pre-registered code never escalated a
  k disagreement at all: its projection multiplied the failing z of k's agreement by a factor above one,
  so a first pass that disagreed on k failed at once, as did one failing the halves' agreement or Rhat.
  The escalated shipped half is written at the move's draw count: its chains continue at five times the
  thinning, and the first pass's draws are thinned by five. The projection reads the first pass the same
  way, its effective draws thinned by five, times five. Before, the projection scaled the unthinned
  series, and the escalated half kept five times the draws: an escalated D1 reference would have needed
  about 45 GB of memory to score.
- **Pooling.** A pooled reference weighs its halves equally: both halves hold the same number of draws
  (above), wherever they could differ the longer is thinned to the shorter, and pooled mass readings
  average the halves' masses rather than summing their counts. Before, an escalated half counted five
  times.
- **S9's fixed-k arms** fix the propensity's k too, at the same quantile of the propensity's k in the
  reference (recorded with every draw, through escalation). Before, the propensity's k was drawn.
- **D2: a reference per pipeline.** The move pipeline's propensity gets its own reference; the move and
  mixed arms are read against it, the shipped arm and the fixed-k arms against the shipped pipeline's. It
  adds about 2.5 core-hours to tranche 1, under the 10 that would have called for one shared propensity
  instead.
- **stan4bart and E3 references escalate.** The smoke's first passes accepted 2 of their 4 (E1 and E3
  probit; E2 and E3 logistic failed), at a tenth of the length; stage 1 ran no references. That is under
  80 percent, so escalation is built: stan4bart reruns the shipped half from the same seed at the
  escalated length and thinning (its first pass, thinned by five, a prefix of the rerun: checked and
  recorded, exact on the forced test), and E3 continues its saved samplers, as the other tasks do. A
  stan4bart rerun holds every iteration's fits, an estimated 8.5 GB for E1 and E2.
- **Acceptance table.** acceptance.rds keeps the per-quantity table: the halves' largest z, joint Rhat,
  tail effective draws and, at a first pass, the shipped half's thinned by five.
- **D1's secondary reading** is built into the reading: the pilot's replications 1 to 10, shipped and move,
  scored as the binary interval study scored them against the designed truth (treatment-contrast coverage,
  M1 to M3, unit and contrast bands), with paired differences; arr1-sep is flagged.
- **Disk.** B keeps no last fit. The audit set is capped at 1 GB of full draws rather than 20 references:
  8, chosen in advance in a seeded order (B1:spambase probit, B2:lalonde logistic, E3 logistic 7 and
  probit 5, S1 probit 6, S5 logistic 6, S6 probit 1, S8 probit 7). Saved samplers drop their control's
  call, which held the data frames: the samplers were 1.7 GB of a D1 reference's 2.7 GB in flight, and an
  S9 sampler falls from 12 MB to 0.2. One reference with 200 MB or more of draws holds them at a time (two
  if the bench box's disk is grown). Projected over tranche 1's schedule: a peak of 6.2 GB and 5.8 GB at
  the end (fits, designs and scores 2.0, compacted references and their readings 2.9, the audit set 1.0),
  against 8.7 GB free once the smoke and stage 1 directories were deleted (their cited outputs copied off
  the box); the queues stop below 2 GB free. With two large references in flight the peak is 7.9 GB, which
  the stop would interrupt. Before, 20 audit references (about 3 GB) and two large references in flight.
- **Memory and workers.** Measured on the box at full length: a D1 reference's fixed-k and scoring jobs
  peak at 14.7 and 14.2 GB, its acceptance and compaction at 7.5 and 8.2, its fits at 5.1; a B1:bank
  reference's jobs at 3.5 or less. Four D1 scoring jobs at once would need 58 GB of the box's 24. Every
  job carries a peak-memory estimate scaled from these by its reference's size (16 GB for D1 and S3
  scoring), and the queue starts the earliest runnable job only while the running jobs' estimates and its
  own stay within 21 GB; a later job passes it only if it leaves the waiting job a worker and its memory
  once the jobs started before it end. Four workers. Alternative: a worker count fixed by the largest job,
  one, about 80 hours of wall time.
- **Cost.** Tranche 1 projects to 83 core-hours from the timing pass's per-sweep times (60 of sampling:
  fits 11, references 37, escalation 10, fixed k 2; plus 15 percent for shared memory bandwidth and 20 for
  scoring and start-up) and 80 from the smoke's job times scaled to full length, against 80 and 77 before
  the amendment (D2's second reference, stan4bart's full reruns, the reference-path control). The
  reference-path control's S6 and S7 cells escalated 2 of 24 references, under the projection's one in
  six. About 30 to 35 hours of wall time at four workers (23 to 27 with two large references in flight).
  This round used about 3 core-hours, 8.5 in all before tranche 1; at the plan's ratio of the whole study
  to tranche 1, the expected total is about 135, under the stop at 211.

## Rulings and pending slots

- **The link in the use tiers**: ruled, dec-B463. A logistic cell counts as the same cell under probit.
- **Full or lean**: the full design, the orchestrator's call under the maintainer's standing permission to
  raise compute ([Considered: a lean variant](#considered-a-lean-variant)).
- **The colleagues' case (D1x)**: asked of the maintainer; pending. If their multiple-treatments problem is
  shared before tranche 1 it runs as a tier 1 reference cell of about 5 core-hours, within the stop; if
  later, as an addition.

## Not verified

Stage 1 settled the first, second and fifth items below ([Pre-registration](#pre-registration)); the
callback is not used, since a callback fit keeps no k.

- That stan4bart's bartcore branch builds against 207f740e's dbarts, that its callback can record the
  scored quantities at every iteration of a long run, and the stan4bart and imputation-loop costs
  (estimated at about 15 and 10 seconds a fit).
- The cost per accepted replication of the rare, weak-signal selections, and that the near-separated band
  can be filled without deterministic draws.
- The per-replication spreads behind [Replications and power](#replications-and-power), and the
  escalation share.
- That references pass the revised acceptance near separation, in at least 80 percent of a cell's
  attempts.
- That the binary interval study's box libraries, fit summaries, references and D1 pilot readings are intact
  on the box.
- The cites to the binary interval study's plan resolve once its records land on bartcore.
