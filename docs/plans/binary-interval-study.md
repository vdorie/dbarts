# binary-interval-study: does the collapsed scale move make binary intervals more accurate at equal running time?

Status: PLANNED 2026-10-17, FOR DECISION. Nothing run beyond the planning probes listed at the end. The
maintainer's calls are under [Before it runs](#before-it-runs); one of them is the open question of
[The multiple-comparisons slot](#the-multiple-comparisons-slot).

agent: one opus implementer writes the harness, runs it on the bench box and writes the finding; one blind
critic of this plan before the run (the k study ran without one and overran its budget 1.5 times).
rng: NEUTRAL. No package code changes; a benchmark harness and two private builds of the move's branch.
window: before 1.0-0; it feeds the collapsed-move plan's open calls (the probit half, a switch) and can
reopen its ruled logistic half.
budget: about 34 core-hours expected, 44 if every diagnosis reference runs, a stop at 50
([Cost](#cost)). Harness about 900 lines of R.

## Summary for the maintainer

You asked for proof that the collapsed scale move helps with what users estimate, not with k: intervals
that are too narrow, and a multiple-comparisons problem, on real binary problems, at equal running time.
What exists is one number: logistic fits on simulated 500-row near-separable data covered the true
probability 79 percent of the time without the move and 99 percent with it. That was at 64,000 draws, not
the default, and on a probit truth that a logistic fit's tails cannot follow, so it says little about a
default fit on real data.

The study compares the shipped sampler with the same sampler plus the move, for probit and for logistic,
on three kinds of problem:

- **Real covariates, made-up outcomes, so the truth is known.** The 22 real binary datasets used in the k
  study, each with its outcome replaced by draws from a probability surface fitted to the real outcome;
  the IHDP infant-development covariates with treatment effects built on them, as in Hill (2011) but with
  a binary outcome; and the ACIC 2016 competition's 58 real covariates with its generated surfaces turned
  into binary outcomes. The last two give individual treatment effects with known truth, including
  designs where every true effect is zero, which is where a multiple-comparisons problem shows.
- **Real outcomes, truth unknown.** The same 22 datasets with their real outcomes, and three real causal
  datasets with binary outcomes (a heart-procedure study with rare deaths, an HIV trial, and the LaLonde
  job-training data). Here each problem gets one very long run as the gold standard, and the question is
  whether default-length runs reproduce its intervals.
- Each fit is read at the default length, at the length the shipped sampler can run in the time the move
  takes (the move costs 14 to 28 percent more per iteration at the default 75 trees), and at a longer
  length, so "just run longer" is compared too.

What is measured: whether 90 and 95 percent intervals contain the truth as often as they claim, for the
predicted probability and for each unit's treatment effect; how wide they are; whether all units'
intervals hold together; how often a unit is declared to have an effect of the wrong sign, or any effect
when there is none; and the error of the posterior mean. The study also asks what goes wrong that the move
cannot touch: whether the exact posterior's own intervals miss (the model or prior), whether the tree
structure mixes too slowly even with the move, and where (rare outcomes, near 0 or 1, units with no
comparable treated units).

The rule, fixed before any result: the move is kept for a link if, at equal running time, it brings
interval coverage at least 2 points closer to nominal on one of three groups (predicted probabilities
with known truth, treatment effects with known truth, real outcomes against the gold standard) and makes
none of them worse by more than 1 point. It is left out if its gain is shown to be under 2 points on all
three, or if it does harm. If neither is shown even after doubling the replications, the probit move is
left out and the logistic move stays as ruled.

Cost: about 34 core-hours, about 10 hours of the bench box with all four cores, less if the halfway
reading already settles it.

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
  the engine's test toggle for the move flipped off by a one-line sed (the patch file kept beside it does
  not apply); on the box, two exports of 207f740e built the same way. The flipped build reproduces the
  base bartcore build bit for bit on seeded probit and logistic fits with an offset, counts and a mask.
  Both builds can therefore be made from the branch for this study, and are its two arms.
- The k study's real-data harness, its 22 datasets and their loaders are in benchmarks/R/binary-hyperprior.R
  and benchmarks/R/uci-binary.R. Coverage practice is that study's (equal-tailed intervals of the per-draw
  probability, coverage distance against 0.90, interior rows split out) and
  [benchmark-surfaces.md](../design/benchmark-surfaces.md)'s (paired margins, a well-mixed reference arm as
  the baseline where the shipped configuration over-covers).

## Problems

### A. Truth known: real covariates, generated outcomes

**A1. The 22 real datasets with a fitted surface as truth.** For each dataset and link, the truth p*(x) is
the posterior mean probability of one fit of that link to the real outcome on the whole dataset (on a
fixed 5,000-row subsample for the five large ones), 4 chains of 2,000 after 500. The truth need not be
converged: it is a fixed surface, whatever chain made it. Each replication draws a new 80/20 split (4,000
and 1,000 rows for the large five, as in the k study) and a new outcome y ~ Bernoulli(p*(x)) on every row;
the fit's held-out intervals are scored against p*. The truth keeps each dataset's character: tictactoe and
mushroom stay nearly deterministic, haberman and birthwt weak. Each link is scored on its own link's truth
(the k study's caveat: a logistic fit on a probit truth measures the two links' tails), with a small
cross-link reading on the side. 22 cells a link, 20 replications.

The truth is a sum of trees, so it favors BART; A2 and A3 have smooth, out-of-class surfaces.

**A2. IHDP with a binary outcome.** The 747 children of Hill's IHDP construction after its deletion rule
(every treated child of a nonwhite mother dropped), 25 standardized covariates, the real treatment. Surfaces
on the latent scale follow Hill's settings: A, a nonlinear control surface exp((x + 0.5) b) against a
linear treated surface x b - w, b sparse; and B, main effects, quadratics and interactions with sparse
coefficients for both arms. Each surface is standardized on the controls and mapped to a probability
p_z(x) = F(a eta_z(x) + c), F the fitted link, c set for a control base rate of 0.3, a = 1 (moderate) or
2.5 (near-separated). The individual effect is p_1(x) - p_0(x). The null cells set p_1 = p_0, so every
true effect is zero. Units whose mother is nonwhite have no treated comparison, so their effects require
extrapolation; coverage is read for them separately. Cells: A moderate, A near-separated, B moderate, null
moderate, null near-separated; 40 replications each, a new surface draw and outcome each replication.

**A3. ACIC 2016 covariates with a binary outcome.** The competition's 4,802 rows and 58 real covariates
(Collaborative Perinatal Project: 3 categorical, 5 binary, 27 counts, 23 continuous), its treatment
assignment and its generated control and treated surfaces, mapped to probabilities as in A2. Settings from
the competition's 77: row 1 (linear surface, high heterogeneity, weak overlap: propensities down to 1e-10),
row 31 (exponential surface, medium heterogeneity, full overlap, low alignment between what drives treatment
and what drives the outcome), row 30 (step surface, high heterogeneity, full overlap), and a null built on
row 31's control surface. At 1,000 units (a new subsample each replication): the three settings and the
null at both signals, base rate 0.3 (8 cells), and row 31 and its null at a base rate of 0.05, moderate
signal (2 rare-outcome cells), 40 replications each. At all 4,802 units: row 31 and its null, near-separated,
20 replications each. 12 cells.

### B. Truth unknown: real outcomes against a gold-standard run

**B1. The 22 real datasets with their real outcomes.** One fixed split per dataset (the k study's first).
The gold standard is a reference run with the move ([Reference runs](#reference-runs)); the arms' 20
default-length runs a dataset differ only in seed, so the spread among them is Monte Carlo error and the
gap to the reference is what poor mixing costs.

**B2. Three real causal datasets with binary outcomes**, individual effects on every unit against a
reference run:
- lindner (twang): 996 coronary-intervention patients, abciximab or not, six-month survival, 26 deaths. A
  rare outcome.
- ACTG 175 (BART): the HIV trial's zidovudine arm against zidovudine with didanosine, 1,054 patients, the
  event indicator (about 27 percent), baseline covariates only. Randomized.
- lalonde (twang): 614 men, job training, earnings in 1978 above zero (77 percent). Observational, weak
  overlap.

### What is offline and what is not

| source | what | how | offline |
|---|---|---|---|
| six R datasets of the k study | pima, biopsy, infert, kyphosis, pbc, birthwt | MASS, rpart, survival, datasets | yes |
| sixteen UCI datasets | adult to sonar | benchmarks/R/uci-binary.R, sha256-checked | cached on the laptop (18 files); the box takes the cache or fetches |
| IHDP | 747 units after deletion | ihdp.RData in the maintainer's npci checkout | yes, local file, copied to the box |
| ACIC 2016 | covariates, 77 settings, generator | the aciccomp2016 package in the maintainer's aciccomp checkout (GitHub vdorie/aciccomp, not CRAN) | yes; installed into the study library in planning |
| lindner, lalonde | real causal, binary outcome | twang | yes |
| ACTG 175 | real trial | BART | yes |
| ACIC 2019 | binary-outcome tracks, high and low dimension | challenge files | no: must be fetched |
| Twins (Louizos et al. 2017) | twin-pair mortality, both potential outcomes observed | CEVAE repository files | no: must be fetched |
| ACIC 2018 / IBM causallib | linked-birth outcomes | challenge files | no: must be fetched, and not confirmed to exist with a binary outcome |

Twins is the one fetchable addition worth having: a real binary outcome where both potential outcomes are
observed for each pair, so individual effects have a real, not generated, truth. It is not in the design.

Duplicate covariate patterns matter for a drawn k: rows that repeat one outcome pull the posterior toward
interpolation (csmq). Measured in planning, the share of rows in repeated covariate patterns with one
outcome is 0.15 in transfusion, 0.12 in spambase, 0.11 in haberman, under 0.03 elsewhere among the UCI
sixteen. Results are split by it.

## Arms

Two builds of the branch tip, each in its own private library:
- **shipped**: the toggle flipped off. Probit runs with its rescaling step, logistic with nothing; this
  is the sampler bartcore ships.
- **move**: the tip as is. Probit runs the step and the collapsed move, logistic the collapsed move.

Package defaults otherwise: 75 trees, chi(1.5, 2), 500 burn-in. Fits that report individual effects use
bartCause's 10 chains, since that is how users of this family of packages estimate them, and fit y on the
covariates and the treatment as bartCause's response model does (no propensity score: that would be a
second BART fit the move also touches); the others use bart()'s 4 chains.

Each fit runs past the default and is read at several lengths, since a longer run's first draws are bit
for bit those of a shorter run with the same seed (checked in planning under both links, and threads 1 and
2 agree). A 4-chain fit runs 500 burn-in and 2,000 draws a chain; a 10-chain fit 500 and 1,000. Readings:

- **default**: the first 500 draws after the 500 burn-in.
- **equal time**: the shipped sampler given the move's running time, from the per-fit times both builds
  record on the same box under the same load. With r the cell's time ratio, the shipped sampler has 1,000 r
  iterations a chain, spent two ways, both read from the same run: as more draws (500 burn-in, 1,000 r -
  500 draws), or as more burn-in (1,000 r - 500 burn-in, 500 draws). The better of the two for the shipped
  sampler, by pooled coverage distance, is its equal-time reading. Choosing after the fact favors the
  shipped sampler, which is the conservative side for the move.
- **long**: all draws of the run (2,000 or 1,000 a chain), for both builds: what running longer buys, and
  how many times the default length the move is worth.

The equal-sweeps comparison (default against default) is reported but does not decide.

## Reference runs

A reference is 16 chains of 5,000 burn-in and 15,000 iterations at thin 15 (16,000 kept draws), with the
move. The move is exact, so the reference targets the same posterior as both arms; the shipped build's
reference is run as well on six B1 datasets and all three B2 datasets, and the two must agree (every
held-out probability's difference within 4.5 combined Monte Carlo standard errors), which also checks the
move's exactness on real data.

Accepted when rank-normalized split-Rhat is at most 1.01 for k and for the mean held-out probability, the
99th percentile of the rows' split-Rhat is at most 1.02, and tail effective draws at the 5th and 95th
percentiles of each row are at least 400. A reference that fails is continued once to twice its length
through its kept sampler; failing again, the problem is reported as having no gold standard reached,
which is itself a finding about tree-level mixing.

References are needed on every B problem (the decision uses them) and, for the diagnosis, on 3 replications
of the A cells that qualify ([Stops](#stops)).

## Measures

All intervals equal-tailed, of the per-draw quantity, at 90 and 95 percent, transformed per draw and
averaged last.

Truth known (A):
- **Pointwise coverage and width** of the probability f(x) at held-out rows (A1), and of each unit's effect
  p_1 - p_0 (A2, A3). Split by rows with true probability inside [0.01, 0.99] and outside, and in A2 by
  units with and without treated comparison.
- **Simultaneous coverage**: one band over all held-out rows or all units, built from the draws (the
  smallest widening of the pointwise intervals that holds 90 or 95 percent of the draws entirely inside,
  Besag et al. 1995); the share of replications whose whole truth vector lies inside. Reported with the
  band's widening factor.
- **Error of the posterior mean**: RMSE against the truth on the probability scale.
- **Sign and significance errors** of individual effects, a unit declared when its 95 (and 90) percent
  interval excludes zero: among declared units the share whose true effect has the other sign; the share
  declared. In null cells every declaration is false: the share of units declared, the share of
  replications with any unit declared (pointwise and by the simultaneous band).
- **Average effect**: coverage of the sample average effect's interval, the estimate most users report.

Truth unknown (B):
- **Reference mass**: the share of the reference's draws inside the arm's interval, per row or unit; 0.90
  for a 90 percent interval when the arm reproduces the posterior. Below it, the arm's intervals are too
  narrow or misplaced. Its floor from Monte Carlo error alone is read by splitting the reference in two.
- **Width ratio** to the reference's interval, and the posterior mean's error in reference standard
  deviations.
- **Simultaneous reference mass**: the share of reference draws entirely inside the arm's band.

Mixing, every fit: split-Rhat and effective draws of k and of the mean held-out probability, wall time.

## What goes wrong beyond the move

The maintainer reports trouble with binary outcomes on real data. Each cause below is read from runs the
design already makes, and goes to the report and TODO whatever the rule decides:

- **The posterior itself misses**: in an A cell with a reference, the reference's own coverage of the truth
  is at least 3 points from nominal. That is the model or the prior, not mixing: nothing about the sampler
  can fix it. Split by near 0 or 1, rare outcome and extrapolation.
- **Mixing the move does not reach**: the move's default reading is at least 3 points of reference mass
  short of the reference, while k has converged. That is the trees, the province of the tree-mixing
  research.
- **Mixing the move does reach**: the shipped reading falls short and the move's does not.
- **The link's tails**: the cross-link reading of A1 (each link fitted to the other link's truth on five
  replications a dataset, shipped build only).
- **Rare outcomes, extrapolation, duplicated rows**: the A3 rare cells and lindner; A2's units without a
  treated comparison and lalonde; the split by duplicate share.

## Decision rule

Fixed before any result is seen. Per link, at equal running time: the move's default reading against the
shipped sampler's equal-time reading.

Three measures of gain, each a coverage distance, the mean over cells of |coverage - 0.90| for 90 percent
intervals:
- **F**: held-out probabilities, truth known (A1's 22 cells);
- **T**: individual effects, truth known (A2 and A3's 11 non-null cells);
- **R**: real outcomes, reference mass in place of coverage (B1's 22 datasets and B2's 3).

Differences are paired by replication or seed; their bounds come from a paired bootstrap over replications
within cells, 2,000 resamples. 95 percent intervals are reported beside and must agree in direction.

**Include the move** for a link when both hold:
1. Gain: F, T or R improves by at least 0.02 (September's bar, the k study's), with the one-sided 95
   percent lower bound on the improvement above zero.
2. No harm: on each of F, T and R the one-sided 95 percent upper bound on any worsening is below 0.01; the
   simultaneous coverage distance (pooled over the A cells) and the average-effect coverage distance are
   not worse by more than 0.02 at the same bound; the pooled RMSE ratio is below 1.02; the null cells'
   share of replications with any unit declared is not higher by more than 0.03 at the same bound.

**Leave it out** when condition 2 fails, or when the upper 95 percent bound on the improvement is below
0.02 on all of F, T and R: its gain in intervals is then shown to be under the bar.

**Otherwise inconclusive**: the replications are doubled once in every cell; still inconclusive, the
numbers go to the maintainer with the recommendation to leave the probit half out (14 to 28 percent of
every default probit fit, not shown to buy interval accuracy) and to keep the logistic half as ruled
(dec-B404, dec-B420: it converges k, and no interval harm was shown).

What each outcome means:
- Probit included: the probit half stays, the collapsed-move plan's large-n options still to be weighed.
  Left out: it is dropped; probit keeps its rescaling step alone.
- Logistic included: as ruled. Left out: dec-B420's consequence, that the move is built, rested on k's
  convergence; under the maintainer's criterion it would not ship, and that goes back to the maintainer.
- Neither: both halves go back, with the diagnosis saying what does make binary intervals wrong.

### Why these counts decide

Planned replications: A1 20 a dataset (440 a link), A2 40 a cell, A3 40 a cell at 1,000 units and 20 at
4,802, B 20 seeds a dataset. With a per-replication spread of the paired coverage difference of 0.05 for F
and R (as the k study's paired standard errors imply) and 0.10 for T (one replication's units share a
surface), the pooled standard errors are about 0.0025 for F, 0.005 for T and 0.002 for R. A true gain of
0.02 is then 4 to 8 standard errors, past the 3.3 needed for 95 percent power at a one-sided 95 percent
bound. A pilot of the first 10 replications of every cell re-estimates the spread; if a measure's projected
standard error passes 0.006, its cells' replications rise to at most twice the plan.

## Stops

- **Before the run** (stage 1): the shipped build differs bit for bit from a build of the current bartcore
  tip on a seeded default fit of either link; or the move build does not differ; or the prefix property
  fails on either build. The first means bartcore has moved binary draws since 3ff28d4c: rebase the branch
  first.
- **Halfway**: when half the replications of every cell and every B problem's reference are in, the rule is
  read at 99 percent one-sided bounds (stricter, to pay for the look). A link that is decided stops; its
  remaining fits are not run.
- **Diagnosis references** run only in A cells where the move's default reading is more than 0.03 from
  nominal coverage, plus four control cells a link chosen now (A1 pima and adult, A2 setting A moderate, A3
  row 1 moderate). Cells already near nominal have nothing to explain.
- **A reference** fails acceptance twice: reported, not extended.
- **Budget**: at 50 core-hours, stop and report what is in.
- **A job** exceeds three times its cell's median time, or a fit returns a non-finite draw: the job is
  logged and skipped; more than 1 percent of jobs skipped stops the run.

## Cost

Per-iteration times on one x86 core, 75 trees, from the k study's long-length fits (probit 0.39 ms at 500
rows, 0.58 at 2,000, 1.19 at 4,000; logistic 0.50, 1.04, 2.06), with prediction at the counterfactual rows
for the effect designs, and the move's cost ratio on its build. A timing probe in planning (747 IHDP units
with their counterfactuals, the default 4 chains of 1,000 iterations, the laptop under load) took 1.7 s
under probit and 2.1 s under logistic, in line.

| stage | what | probit | logistic |
|---|---|---|---|
| 1 | builds, checks, timing pass on the idle box | 0.5 | 0.5 |
| 2a | A1 truths and fits (880), cross-link reading (110) | 1.8 | 2.7 |
| 2b | A2 fits (400, 10 chains) | 0.9 | 1.3 |
| 2c | A3 fits (960, 10 chains) | 2.7 | 4.9 |
| 3a | B1 fits (880), B2 fits (120) | 1.8 | 2.8 |
| 3b | B references: 25 with the move, 9 shipped | 1.9 | 2.8 |
| 4 | diagnosis references, A cells that qualify (up to 117 a link) | up to 7.6 | up to 11.9 |
| | total, all diagnosis references | 17.2 | 26.9 |
| | total, expected (half of them) | 13.4 | 21.0 |

Core-hours. Expected about 34, about 10 hours of the box at one job a core (four jobs share memory
bandwidth: 15 percent added); at most 44 with every diagnosis reference; with the halfway stop firing on
both links, about 25. The stop is at 50.

Memory: the largest jobs (A3 at 4,802 units, 10 chains of 1,000 draws, fits at both treatment values) hold
about 0.8 GB of draws and twice that at peak; they run at most two at a time. Disk: individual-effect draws
kept as single precision, thinned to 500, at the default and equal-time readings of every effect fit, about
12 GB, so an estimand added later needs no new fits ([The multiple-comparisons slot](#the-multiple-comparisons-slot)).
Per-fit summaries are small.

## Scripts to write

- **benchmarks/R/binary-interval-study.R**, tracked, about 900 lines: the problem registry (A1 to A3, B1
  and B2, their cells and seeds), the generators (A1's truths, Hill's two surfaces on IHDP, ACIC 2016
  through aciccomp2016's generator, the link mapping, the null and rare variants), the fit runner with its
  readings, every measure above, the reference runner with its acceptance and one continuation, and actions
  `jobs`, `truth`, `fit`, `reference`, `pilot`, `interim` and `summary`, the last two computing the rule
  exactly as written here. It sources benchmarks/R/uci-binary.R; the six R-dataset loaders move from
  binary-hyperprior.R into a file both scripts source, order and seeds unchanged. The IHDP file's path
  comes from an environment variable.
- The box driver, in scratch: a queue of four workers, each job running both builds in alternating order,
  jobs in random cell order so the halfway reading covers every cell; the largest jobs limited to two at a
  time.
- The flipped build: a correct one-line patch setting the move's toggle default to false, applied to an
  export of the tip (the patch kept in the move's scratch does not apply).

Outputs stay in scratch/bis; the finding goes to docs/design/binary-interval-study.md.

## Order of work

1. Write the harness on the laptop; smoke it with one replication and 50 iterations a cell, both links,
   both builds.
2. On the box: both builds from 207f740e, `--preclean`, aciccomp2016, the data. Stage 1's checks and the
   timing pass with the box otherwise idle.
3. The pilot: 10 replications of every cell and every B problem's reference; sizes reread.
4. Stages 2 and 3 to the halfway point; the halfway reading.
5. The rest of stages 2 and 3 for an undecided link; stage 4 for the qualifying cells.
6. The summary, the finding, a ledger entry for the calls made running it, and the result to the maintainer
   with the collapsed-move plan's open calls.

If the move lands on bartcore before this runs, its two builds are bartcore's tip and the tip flipped, and
stage 1's first check is moot.

## The multiple-comparisons slot

**OPEN, awaiting the maintainer**: what multiple-comparisons problem have you seen with binary fits?

How the design takes each likely answer:
- Many units' effects declared nonzero when there are no effects, or when they are homogeneous: the null
  cells measure it directly (share of units declared, share of fits with any unit declared). Nothing to add.
- Intervals read jointly, all units at once: simultaneous coverage. Nothing to add.
- Subgroups searched for effects after looking (by covariate level or quantile): subgroup average effects
  and their joint coverage, computed from the stored effect draws. No new fits.
- Pointwise coverage below nominal across many replications or units: the pointwise measures. Nothing to
  add.
- Many outcomes or many models on the same data: a new design, several outcomes on one covariate set. About
  2 core-hours a link for a five-outcome cell at 1,000 units.
- Something else (a multilevel fit through stan4bart, bartCause's propensity-adjusted fits, variable
  importance): a new arm or design, brought back with its cost before the run.

The run waits for the answer only if it names a new design; otherwise the answer is absorbed into the
summary.

## Before it runs

The maintainer's calls, each with a recommendation:

1. **The multiple-comparisons question**, above.
2. **The bar.** Recommended: September's 2 points of coverage distance for a gain, 1 point for harm.
   Coverage distance treats intervals too wide like intervals too narrow; the alternative counts only the
   shortfall below nominal, which matches "too narrow" but would credit a move that over-covers. The
   logistic number so far was over-coverage (0.99 against 0.90).
3. **How the shipped sampler spends equal time.** Recommended: the better of more draws and more burn-in,
   which favors the shipped sampler. The alternative is the one a user would most likely choose, more
   draws, fixed now.
4. **The A1 truth.** Recommended: a fitted BART surface, which keeps each real dataset's character but
   favors BART, with A2 and A3 as the out-of-class check. The alternative adds a smooth fitted truth (an
   additive model with interactions) on the 22 datasets: about 4 core-hours more.
5. **Fetching data.** Recommended: no fetch beyond the UCI cache. Twins would add a real binary outcome with
   both potential outcomes observed, about 3 core-hours, if the maintainer allows the download.
6. **Chains for effect fits.** Recommended: bartCause's 10, as users of this family of packages fit them.
   The alternative, bart()'s 4, saves about 2.5 core-hours.
7. **The budget**: about 34 core-hours, about 10 hours of the box, stop at 50.

## Probes run in planning

On the laptop under a load of 15 to 50, private library scratch/libs/bis (the branch tip and aciccomp2016),
one job at a time:
- the build of 207f740e installs; aciccomp2016 and aciccomp2017 install offline from the local checkout;
- data loads: ACIC 2016 (4,802 by 58, 858 treated in its first setting's first replication, the generator
  0.25 s), IHDP (747 after the deletion, 139 treated), lindner (996, 26 deaths), lalonde (614, 77 percent
  with earnings), ACTG 175 (2,139, 4 arms), the 16 UCI datasets from the cache;
- one default fit under each link on IHDP with a binary outcome (1.7 s and 2.1 s);
- a 2,000-draw run's first 500 draws equal a 500-draw run's, bit for bit, under both links, and threads 1
  and 2 draw alike;
- duplicate covariate patterns in the UCI datasets.

## Not verified

- Whether bartcore has moved plain binary draws since the branch's base (stage 1 checks it).
- The x86 time of the effect designs and of 10-chain fits; the cost table scales the k study's timings.
- The per-replication spread the counts are sized on (the pilot rereads it).
- That every reference reaches acceptance at the planned length, probit near separation least of all.
- The memory peak of the largest jobs on the box.
