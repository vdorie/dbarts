# nbinom-default-study: is there a better default prior for negative binomial fits?

Status: PLANNED 2026-10-29 under dec-B432, revised the same day from its blind critique, and amended
after the critique's re-check (verdict: build stages 0 and 1;
[Amendments before results](#amendments-before-results)). Stages 0 and 1 ran the same day
([Stages 0 and 1](#stages-0-and-1)): the move build passed its checks, and as written no candidate passed
Part A. The maintainer ruled the criterion (dec-B502): five candidates go to the fits, under the
[Amendment after stage 1](#amendment-after-stage-1). Stage 2 ran the same day ([Stage 2](#stage-2)): the
long reading's check passed, and the main run re-prices at about 60 core-hours before confirmation. The
orchestrator ruled on its report the same day ([Amendment after stage 2](#amendment-after-stage-2)). C4's
and C5's scratch build passed its drift tests, Part B is in the harness, revised on the Part B review's
findings, and the smoke ran clean on its own seeds ([Box readiness](#box-readiness)): about 48 core-hours
expected before confirmation, the stop at 110. Pre-registered the same day, before any study fit was
scored ([Pre-registration](#pre-registration)). Runs after the merge to main and before 1.0-0.

agent: the blind critique (opus) ran; this revision answers it. One opus implementer writes the harness,
runs it and writes the finding; an opus review of the harness, the scratch build and the smoke at stage 4,
since the harness carries the statistics; one opus review of the main reading before it goes to the
maintainer.
rng: NEUTRAL. No package code changes: a benchmark harness and its outputs. The long reading runs on a
scratch build carrying the collapsed scale move's negative binomial half (an evaluator and the family in the
move's predicate, about 15 lines on branch collapsed-move-main), never landed. A finalist the package cannot
express today is checked at the default length on a scratch build carrying only its k draw.
window: after the merge to main, before 1.0-0 (dec-B432). A change of default it supports goes back to the
maintainer; building it is a separate, posterior-changing slice.
budget: harness about 1,100 lines of R (benchmarks/R/nbinom-default-study.R; stages 0 to 2 and Part B,
written, are about 5,300 as formatted), reusing the first study's prior predictive draw. Compute about 48
core-hours expected before confirmation, stop at 110 ([Box readiness](#box-readiness)).

## Summary for the maintainer

You kept the negative binomial default, k drawn under chi(1.5, 2), judged the study behind it flawed and the
default improvable (dec-B431), and ruled that a redesigned study runs before 1.0-0 with its criterion fixed
first and in absolute terms, count data with wide spreads, chains long enough for k to mix, and candidates
beyond the prior on k (dec-B432). This is that study.

**What it can show.** On data of the sizes here the data, not the prior, set k: its posterior sd on the log
scale is 0.10 to 0.20 on real data, against 0.35 to 0.80 under the candidate priors. So candidates that
differ in the prior's tails will mostly give the same intervals. What a change of default buys is a prior
whose datasets look like count data, not better intervals; on fits, the study checks that a candidate loses
nothing. Its compute goes where candidates can differ: small samples with wide spreads, and narrow or null
data, where a prior with a thinner tail in k widens intervals the default keeps tight (about 20 percent
wider on null data in the critique's probe).

**The criterion, ruled by you after stage 1.** Two parts.

- Part A, the prior alone, in absolute terms. Recommended: in at most 1 prior dataset in 100 should the rows'
  mean counts span more than a factor of about 10^8 (a spread of the log mean across rows, S, above 4); at
  most 1 in 100 should hold a count the package refuses (above 10^6); and each real dataset in the study
  should be no rarer than 1 in 100 under the prior, measured over the rows whose mean the data identify.
  Stage 1 computes Part A for every candidate (about one core-hour, no fits). You then rule on the criterion
  with that table in hand, before any fit is scored; the fits stay blind until then.
- Part B, what fits must not lose: interval coverage of what users estimate from count fits (the mean count
  at a row, mean counts by group, rate ratios of a binary covariate, a new count), on truths with known
  values built on real covariates, read as distance from nominal on either side, with the house's bar of
  0.02, and a guard on interval width on narrow and null data. That guard, 10 percent, may exclude the
  thinner chi priors: on null data the critique's probe found their intervals about 20 percent wider than
  the default's.

A candidate that passes Part A replaces a default that fails it if it loses nothing on Part B; were the
default to pass, a candidate would have to be better by the bar.

**Compared**: the default; k fixed at 1.5 and at 2; k drawn under chi with 3, 4, 5, 7 or 10 degrees of
freedom at the default's median; the default's form at two-thirds the width; the default with its tail cut
at three places (k kept above a floor, so the forest's spread has a ceiling); a half-normal prior on the
forest's spread at two widths; a spread scaled to the data, as gaussian's sigest is; and a shape prior that
puts more weight on small shapes.

**Chains long enough for k to mix.** Under today's sampler k's autocorrelation time was 440 to 8,000 sweeps
on the eight datasets probed. The collapsed scale move, built for the binary links and not yet for
negative binomial, was given its negative binomial half on a scratch build: it samples the same posterior
(it agreed with today's sampler wherever today's converged) and cut k's autocorrelation 5 to 17 times on
narrow and moderate data and 2.5 to 5 times on wide. The study's long reading is long runs on that build,
one per replication, each reweighted to every candidate prior on the spread: exact, paired across
candidates, and about a third cheaper than integrating k over a fixed-k grid once the grid's error is
counted honestly. On two datasets the trees themselves mix slowly; there runs are capped and flagged. What
users get from a default-length run is read on today's sampler, for the default and the chosen candidate.

**Data**: twelve real count datasets with spreads from about 0.15 to about 2.4 over identifiable rows
(named below), with truths fitted to each by a fixed rule and widened to spreads of 2 and 3; small samples
of 30 to 100 rows from real covariates at spreads of 2 and 3; designed cells for exposure offsets, mostly
zeros, high counts, no signal, narrow spread and 5,000 rows; and the real outcomes, scored on held-out
counts. Cells are weighted equally by spread band: narrow, middle, wide, and small-sample wide.

**Stages**: harness and checks; the prior check of every candidate; your ruling on the criterion; the
check of the long reading; a smoke; the main run; the default-length check of the chosen candidate; a
confirmation only if the call is close; the finding to you. About 100 to 150 core-hours expected, stop at
160 (amended: the candidates likely to pass Part A are the thin ones, which need their own runs more
often). The laptop runs everything up to the smoke, and the main run's narrow, middle and
small-sample cells, on two cores; the box runs the smoke and the wide cells after the k study's first
tranche.

## Goal

After this study each candidate default for `family = "nbinom"` carries a pass or fail on each part of a
criterion the maintainer fixed before any fit was scored, with the numbers behind it, and either one
candidate is recommended to the maintainer as the new default or the default is recommended kept, with what
keeping it costs.

## What is known

- **The prior as built.** A row's log mean is c + f(x) + o: c is the data's log rate, a transform the fit
  sets (the log of the total count over the total exposure), o the offset, f the forest
  ([nbinom-log-mean.md](nbinom-log-mean.md), [Leaf prior](nbinom-log-mean.md#leaf-prior)). Given k the
  forest's sd at every x is s = A / k, A = 3 ([`defaultLeafScale`](../../R/model.R)), whatever the tree
  count. k is drawn under chi(1.5, 2) ([`resolveLeafHyperprior`](../../R/model.R),
  [`drawsLeafKByDefault`](../../R/spec.R), the conjugate draw [`ChiKHyperprior`](../../src/bartcore/model.hpp)),
  the same prior as sd = invchi(1.5, 1.5): median k 1.91, median s 1.57, P(k > 4) 0.085. The shape r is
  drawn on the 13-point grid 1 to 50 under gamma(2, 0.1) weights
  ([`NBShapePrior`](../../src/bartcore/model.hpp)): 15 percent of its mass at r of 3 or less, 46 percent at
  10 or more.
- **The refused count.** A count above 10^6 is refused at creation and at every response swap
  ([`maximumCount`](../../src/bartcore/model.hpp), dec-B202): the exact Polya-Gamma draw costs one draw per
  trial, about 160 ns, every row and sweep, and the shape's count histogram is sized by the largest count.
- **The first study** ([nbinom-prior-predictive.md](../design/nbinom-prior-predictive.md)). Under the default
  S has median 0.93, 90th percentile 3.1, 99th 15; at a mean count of 10, 5 to 10 percent of prior datasets
  hold a refused count. No candidate fitted better beyond fold noise; k fixed at 2 covered 0.65 against 0.92
  on wide data at n = 30, its one clear signal. Its flaws are listed in
  [Limits of this study](../design/nbinom-prior-predictive.md#limits-of-this-study) and each is answered
  here: the criterion is absolute and ruled before any fit is scored ([The criterion](#the-criterion)); the
  cap is the package's, 10^6; the candidates and the rule are fixed before any fit
  ([Candidates](#candidates), [The decision rule](#the-decision-rule)); the data are wide and include offsets
  and correlated columns ([Data](#data)); k mixes in the long reading
  ([The long reading](#the-long-reading)); and the shape's prior and a data-scaled spread are among the
  candidates.
- **k's mixing.** On datasets drawn from the default prior at n = 150, k's autocorrelation time has median
  191 sweeps and 99th percentile 3,253 ([k-mixing-pg-families.md](../design/k-mixing-pg-families.md)). On
  the probed data it is 1,140 to 8,000 sweeps where the spread is 1 or more and 440 to 620 on narrow data,
  where k's posterior lies mostly above 4 (Traffic, Insurance, null data: 73 to 98 percent)
  ([Probes](#probes)).
- **The posterior is the data's.** The posterior sd of log k is 0.10 to 0.20 on the probed real data; the
  prior sd of log k is 0.80 under chi(1.5, 2), 0.40 under chi(4) and 0.35 under chi(5). Where the data say
  little (null data, small samples) the prior shows: on null data at n = 500 the critique's probe found
  the default's median k at 6.2 and chi(4, 1.03)'s and chi(5, 0.9)'s at 4.1 and 3.8, with every prior
  covering 0.99 and the thinner tails' intervals about 20 percent wider. It is the only fit-based evidence
  seen before the criterion is ruled, disclosed so that the ruling is made knowing it; it is why Part B has a
  width guard.
- **Cost grows with the counts.** A sweep of today's sampler costs about 0.16 microseconds per unit of y + r
  summed over rows, r the drawn shape, plus about 0.9 microseconds per row, at 75 trees; the scratch move
  build costs 17 to 24 percent more a sweep ([Probes](#probes)). A wide cell's cost is set by its count
  total.

## The criterion

Two parts: A, the prior predictive, decided on prior draws alone in stage 1; B, the posterior, decided on
fits. Both are the maintainer's judgments. This section states the recommendation; the maintainer rules on
it after stage 1 with stage 1's full Part A table, before any Part B fit
([Open calls for the maintainer](#open-calls-for-the-maintainer), 1). The ruling is quoted in the
pre-registration. The finding applies it and nothing else: other criteria appear there only as a sensitivity
section after the recommendation, labeled as such.

### Part A: what a prior dataset should look like

**Convention.** A prior dataset is drawn at a cell's design with the center set as the fit would set it:
c is chosen so that the dataset's expected total count is the cell's mean count times its total exposure,
as the data-derived c of a fit to that dataset makes it. So a cell's label is the dataset's realized mean
count, whatever its spread.

Per prior dataset: S, the sd across rows of f(x) (the offset excluded: it is data, not prior); S_id, the same
over the rows whose expected count (exposure included) is at least 0.01; and the largest count. Cells: each
real dataset's own design (its covariates, n, offset and mean count), and a simulated grid of n 50, 500,
2,000 and 10,000, mean counts 0.1, 1, 10 and 100, and three designs (5 correlated normal columns; 20
uniform columns; 5 normal and 5 binary columns with an exposure offset spread over three orders of
magnitude). 20,000 prior datasets a cell. A candidate passes A when, in every cell:

- **A1, spread.** P(S > 4) is at most 0.01. S = 4 means the middle 98 percent of a dataset's rows have mean
  counts spanning a factor of exp(2 x 2.33 x 4), about 10^8; the widest real data here span about 10^5 over
  their identifiable rows. A prior that goes past 10^8 in more than 1 dataset in 100 spends its mass where
  no count data have been seen.
- **A2, fittable.** P(a count above 10^6) is at most 0.01. A default should not generate, in more than 1
  dataset in 100, data the package refuses. Under the convention above the largest mean count is at most
  the cell's total count, so A2 can bind only in cells whose total reaches 10^5 or more (n 2,000 at a mean of
  100, and n 10,000 at 10 and 100); it does not read A1's tail a second time. Every candidate is scored on
  it in stage 1.
- **A3, room.** For each real dataset, the top of its range of S_id across four fits (a quadratic NB GLM,
  the additive NB GAM of the truths' rule, and, where it has two numeric covariates, s() and te() smooths of
  them; [Real datasets](#real-datasets)) is at or below the prior's 99th percentile of S_id at its own
  design and mean count. S_id leaves out rows whose mean only the prior sets: a fitted mean of 10^-5 in
  an all-zero corner is not something the data ask the prior to reach. Over those rows the spread is
  robust to the fit: within 7 percent of its midpoint across the four fits on every dataset, where over
  all rows crimtab's ran from 2.75 to 3.80. The level, 1 in 100, is A1's: the twelve datasets were chosen
  for being wide, and the widest of a dozen so chosen should still not be a 1-in-100 dataset under a
  default meant for count regression in general. With A1, it places the prior's 1-in-100 upper tail
  between the widest real data and a span of 10^8.

Reported beside, not read: the quantiles of S and S_id, the share of datasets with more than 90 percent
zeros, the count total (cost), each real dataset's percentile, and the narrowest real dataset's (Traffic)
percentile at the low end, which A does not bound and the width guard of Part B answers instead.

By the first study and the probes the default fails A1 (5 to 8 percent above 4 at the probes' six designs).
A2 and A3 are known for no candidate yet: stage 1 computes them.

### Part B: what fits must not lose

Estimand classes, what users report from a count regression:

- **F**, the log mean at each row, training and held-out (equivalently the mean count); pointwise.
- **G**, the mean count averaged over a covariate-defined group (a factor's levels, a numeric covariate's
  quartiles); one scalar per group.
- **R**, rate ratios of a binary covariate where the design has one: the per-row log rate ratio
  f(x, 1) - f(x, 0) and the sample-average rate ratio, the ratio of average mean counts at 1 and at 0.
- **P**, a new count at a held-out row: 90 percent predictive coverage by randomized PIT,
  u = F(y - 1) + V (F(y) - F(y - 1)).
- Rows whose true expected count is below 0.01 are scored apart, in F and P, and reported: their coverage is
  decided by how far each prior lets the log mean fall, not by the data (dec-B462's lesson from the binary
  pilot). The shape r is reported, not read.

Measures, per cell, class and candidate: coverage of 90 percent intervals (50, 80 and 95 read beside); the
cell's gap, its distance from nominal on either side corrected for noise as
[Decision rule](binary-interval-study.md#decision-rule) defines it (dec-B455); the interval score and mean
width; RMSE of the log mean; and on real outcomes the held-out log predictive density per row.

**Bands and weights.** Each cell belongs to one band by its truth's S_id and n: narrow (S_id below 0.5, and
the null cell), middle (0.5 to 1.5), wide (1.5 or more, n above 130), small-sample wide (n 130 or less, S_id
2 or more). A class average is the mean of the four bands' averages, cells equal within a band.

**The long reading** is each candidate's posterior from runs in which k mixes
([The long reading](#the-long-reading)): the prior's own effect, as dec-B432 asks. **The default reading**
is a run of today's sampler at bart()'s defaults (4 chains of 500 burn-in and 500 kept, 75 trees), what
users get. A candidate's draws can move the drawn shape and so its cost; each candidate's time is reported
from its own runs.

Bars, the house's (September's, reused by the k study,
[The bars](../design/binary-k-fixed-or-drawn.md#the-bars)): 0.02 of coverage gap; 0.01 nats per row of
held-out log score. Bounds: an average is "worse by the bar" when the upper end of its 95 percent interval
(bootstrap over whole replications, paired) reaches the bar; a single cell when its point estimate exceeds
the bar by 1.96 of its standard errors.

### The decision rule

1. **A**, from stage 1, under the criterion the maintainer ruled: each candidate passes or fails. A
   candidate that fails A is reported on B where it costs nothing, never recommended.
2. **Non-inferiority, long reading**: in each of F, G, R and P, the candidate's band-weighted average gap
   less the default's is below 0.02 at its upper bound; no cell is worse by the bar; the held-out log score
   on real outcomes is not worse by 0.01 nats per row at its upper bound; the interval score is not worse by
   more than 5 percent in any class; and the **width guard**: on each narrow-band cell, F's interval score is
   not worse than the default's by more than 10 percent at its upper bound. Cells whose truth has a shape
   below 1 are a guard of their own: none worse by the bar.
3. **If the default fails A** (expected): among candidates that pass A and step 2, the one with the smallest
   band-weighted average gap over F and R is recommended. Candidates within 0.005 of it count as tied, and a
   tie goes first to a candidate that passes A with margin (A1's and A2's probabilities at most half their
   bounds, A3 at the prior's 98th percentile), then to the smaller code change ([Candidates](#candidates)),
   then to the smaller interval score over the narrow band. If none passes, the default is recommended kept,
   with what passing A costs in coverage.
4. **If the default passes A**: a candidate is recommended only if it also improves F or R by at least
   0.02 at its lower bound; otherwise the default is kept.
5. **The default reading guard**: the recommended candidate's default reading is not worse than the
   default's by the bar in any class (step 2's bounds). If it is, the next candidate by step 3's order is
   checked. A candidate better than the default only at the default reading is reported as a mixing
   finding for TODO k-mixing-pg-families, not as a better prior.
6. **Close calls are confirmed**: if the recommended candidate's step 2 margin in any class, or its lead
   in step 3 or 4, is under half the bar, a confirmation on fresh replications (a third of the main run's
   count, every cell) repeats steps 2 to 5 on the long and default readings, and decides.

### Stored outputs

Per prior dataset: S and S_id at unit spread and the 20 largest values of f, from which S, S_id and the
largest count under any prior on the spread follow by rescaling (S is proportional to the spread; S_id and
the largest count are recomputed from stored full draws at the real designs). Per replication and run: per
draw k, r and the run's log weights for every candidate; per row the estimands' draws thinned to 500, as
single precision, so that any prior on the spread or the shape can be applied later by reweighting; G and R
draws; P's cdf at y and y - 1 per draw; the truth; effective draws and split-Rhat. About 6 GB.

## Data

### Real datasets

All available offline: base R's datasets and recommended packages, and CRAN packages installed on the
laptop, exported once to .rds files for the box. S over all rows and S_id over rows with a fitted expected
count of at least 0.01, offset excluded, from four NB fits by the rule of
[Truths on real covariates](#truths-on-real-covariates); theta is the additive GAM's shape
([Probes](#probes)).

| dataset (package) | n | total count | zeros | S_id across fits | S, all rows | theta | note |
|---|---|---|---|---|---|---|---|
| roaches (rstanarm) | 262 | 6,720 | 0.36 | 1.20 to 1.32 | same | 0.33 | exposure offset |
| solder (survival) | 900 | 4,977 | 0.32 | 1.36 to 1.37 | same | 4.2 | factors only |
| aids (boot, dud == 0) | 465 | 6,215 | 0.21 | 1.63 to 1.71 | same | 5.9 | |
| mite, long form (vegan) | 2,450 | 9,800 | 0.57 | 1.55 to 1.58 | same | 0.44 | 35 species by 70 cores |
| Arabidopsis (lme4) | 625 | 18,727 | 0.20 | 1.07 | same | 0.60 | factors only |
| grouseticks (lme4) | 403 | 2,567 | 0.31 | 1.36 to 1.37 | same | 0.77 | |
| epil (MASS) | 236 | 1,950 | 0.10 | 0.75 to 0.77 | same | 3.2 | |
| Insurance (MASS) | 64 | 3,151 | 0.02 | 0.31 | same | above 10^5 | holders offset; Poisson |
| crimtab, long form (datasets) | 924 | 3,000 | 0.67 | 2.12 to 2.43 | 2.75 to 3.80 | 0.75 | finger length by height; 11 to 33 percent of rows below 0.01 |
| BCI, long form (vegan) | 11,250 | 21,457 | 0.60 | 1.79 | same | 1.5 | 225 species by 50 plots |
| Traffic (MASS) | 184 | 3,965 | 0.00 | 0.15 to 0.16 | same | 11 | |
| quakes (datasets) | 1,000 | 33,418 | 0.00 | 0.49 | same | 20 | stations reporting |

Five have a shape below 1, the shape grid's floor (TODO approximate-distributions). Excluded from fits and
reported against each candidate's percentiles: Migration (carData; S_id 1.60; its counts total 830,000,
about 1.5 hours a replication).

### Truths on real covariates

For each real dataset, two truths, fixed once from the real outcome and stored:

- **Additive**, by a rule fixed now: an NB GAM (mgcv's nb() family, REML), with a smooth s(x, k = min(10,
  distinct values - 1)) of each numeric covariate with more than 10 distinct values, every other covariate as
  a factor, the offset as given, and no tensor, bivariate or by-terms. Its spreads are the table's
  "additive" fits ([Probes](#probes)).
- **BART-fitted**: the posterior mean log mean of a long dbarts fit at k fixed at 1 (s = 3, wider than every
  candidate's center, so the truth is not one candidate's own shape). It carries interactions the additive
  rule leaves out; where the fitted shape sits at its floor it also carries dispersion as spread, which is
  reported beside its S_id.

Counts are drawn as rnbinom(size = theta, mu = exp(truth + offset)) with theta the additive GAM's shape,
floored at 1 so that the model can reach the truth; for the five datasets whose shape is below 1, an extra
as-fitted cell at the GAM's own shape is read as the shape-below-1 guard. Held-out rows are a fixed fifth of
each dataset's rows.

Spreads: **as fitted** (both truths), and the additive truth **widened** to S_id of 2 and of 3 where it is
narrower, by scaling the truth's departure from its mean. A widened truth whose expected count total
exceeds 20,000 or whose expected largest count exceeds 10^5 has its center lowered until it does not, and
carries the resulting mean count as a label.

### Small-sample wide cells

Where the prior reaches the posterior. Rows subsampled afresh in each replication from the covariates of
roaches, aids, epil and quakes, at n of 30 and 100, with the additive truth widened to S_id of 2 and 3 at a
mean count of 5 and shape 3: 16 cells. Insurance's widened truths (n 64) and D4 at S of 2 and 3 join the
band.

### Designed cells

Simulated covariates, each isolating one feature; truths are smooth functions with interactions, scaled to
the stated S:

| cell | n | feature | S | mean count | shape |
|---|---|---|---|---|---|
| D1 offset | 500 | exposure log-uniform over 1 to 1,000 | 1 | rate 0.05 | 3 |
| D2 zeros | 1,000 | mostly zeros | 1 | 0.2 | 2 |
| D3 high | 300 | high counts | 1 | 100 | 10 |
| D4 small | 60 | small sample | 2 and 3 | 5 | 3 |
| D5 null | 500 | no signal | 0 | 5 | 5 |
| D6 large | 5,000 | many rows | 1.5 | 2 | 3 |
| D7 effect | 1,000 | binary treatment, heterogeneous log rate ratio 0.3 to 1 | 1 | 5 | 3 |
| D8 narrow | 1,000 | narrow spread | 0.3 | 5 | 3 |

Every designed cell has 5 correlated normal covariates (correlation 0.5 between neighbors) and one binary
column; D7's treatment depends on them, as an observational design would.

### Real outcomes

Each real dataset's own outcome, scored by 5-fold held-out log predictive density per row and P's coverage
of the held-out counts, two repetitions. No truth is involved, so it enters only P and the log-score guard.

## Candidates

Every candidate below but C6 is a prior on the forest's spread s = 3 / k, so each is read from the same
long runs by reweighting ([The long reading](#the-long-reading)), and a further one costs almost nothing to
score; each form is sampled at several widths. A different leaf scale A is not a candidate: the posterior
depends on A only through A / k, but as a change of default a new A would also move every fixed k and chi()
a user states, so it is not an implementation of C3. The tree prior and tree count are not candidates: the
tail is in s, and S's dependence on the tree prior is reported from stage 1.

| | candidate | s's prior | why | code change if chosen |
|---|---|---|---|---|
| C0 | default | k ~ chi(1.5, 2), s = 1.5 / chi(1.5), median 1.57 | | none |
| C1 | k = 1.5, k = 2 | s fixed at 2 or 1.5 | no tail at all, no k to mix | a family default for k in place of [`drawsLeafKByDefault`](../../R/spec.R)'s flag; R about 30 lines, help, prior-defaults.md, tinytest; fits lose their k channel |
| C2 | chi(3, 1.25), chi(4, 1.03), chi(5, 0.9), chi(7, 0.75), chi(10, 0.62) | the default's median k, 1.9, thinner tails both ways | the first study's best-looking direction, now chosen before fits | as C1 |
| C3 | chi(1.5, 3) | the default's form at two-thirds the width | narrows bulk and tail together | as C1 |
| C4 | chi(1.5, 2) with k at least 0.6, 0.75 or 1 | the default, s capped at 5, 4 or 3 | cuts only the tail A1 objects to; leaves the narrow side exactly | [`ChiKHyperprior`](../../src/bartcore/model.hpp) gains a floor, k^2 drawn from its truncated gamma by inversion (exact, constant cost; rejection would loop where the conditional sits at the floor, as on wide data it does); [`chi`](../../R/model.R) and invchi gain a bound, carried by the bridge and the internal parameters; no ABI event (dbarts.h carries no hyperprior); about 150 lines and an exact gate |
| C5 | half-normal on s, median 1.3 or 1.57 | thin tail without a hard cap; density largest at s = 0 | the standard weakly informative form | a non-conjugate k step (slice on log k), a constructor, the bridge, the internal parameters, docs that say invchi is the only form; no ABI event; about 300 lines and an exact gate |
| C6 | the shape's grid weights from a continuous log-uniform on [1, 50] | weight of each grid point the log-width of its cell, cells split at the geometric means of neighbors: 0.0886, 0.1404, 0.0886, 0.0653, 0.0518, 0.0601, 0.0653, 0.0518, 0.0518, 0.0653, 0.0886, 0.1171, 0.0653 for r = 1 to 50; 32 percent at r of 3 or less (15 under gamma(2, 0.1)) | real counts mostly have small shapes ([Real datasets](#real-datasets)); crossed with the default and with step 3's choice | the weights in [`NBShapePrior`](../../src/bartcore/model.hpp) and its mirrors in benchmarks/R; help; about 40 lines |
| C7 | data-scaled | k ~ chi(4, c), c setting the median of s to 1.7 times the spread of a Poisson GLM's linear predictor (main effects, offset in the fit, out of the spread), the spread floored at 0.25 | the sigest precedent: gaussian's default reads its scale from an lm fit; 1.7 makes the prior's median S about the GLM's spread | the GLM fit in the default's resolution, as sigest's lm; help; tests; R about 80 lines |

C6 is scored by per-draw weights on r computed in the harness at run time, against the weights above, fixed
now. C7's Part A depends on the data: it is read at each real design at that dataset's own GLM spread, and
reported on the simulated cells at spreads of 0.25, 1 and 2.

Every change of default is posterior-changing for nbinom fits: the exact gates' nbinom entries, the nbinom
equivalence scenario and SBC arm re-derived, and a design note, per the
[RNG classes and their gates](README.md#rng-classes-and-their-gates).

Not candidates, with the reason: a real or sub-1 shape (not built, TODO approximate-distributions); the
grid's cap at 50 (Poisson-like data at large means; a separate question); the center c (data-derived; moving
it is a different model, not a prior).

## The long reading

**What it is.** Per replication, one long run of 4 chains on a scratch build carrying the collapsed scale
move's negative binomial half, under C0's prior. The move draws k jointly with a rescaling of every leaf, the
latents integrated out (the k study's move arm, docs/plans/collapsed-scale-move.md on branch
collapsed-move-main); its negative binomial half, which that plan leaves as a door (its open call 4), is an
evaluator of the log-likelihood at scaled fits and the family in the move's predicate. It is exact: it
changes how the chain moves, not the posterior.

**Every spread prior from one run.** Only k's prior differs between C0 and C2 to C5 and C7, so each one's
posterior is C0's reweighted by pi_c(k) / pi_0(k), per draw, exactly; C4 is the truncation. Coverage is read
exactly from the weighted share of draws below the truth. Where a candidate's weights keep fewer than 30
percent of the run's effective draws (expected on narrow, null and small-sample cells, where the prior shows,
and for C4 where the posterior sits at its floor), that candidate gets its own run on that replication; C4's
own run, and C5's where needed, use a scratch build with its k draw. C1 gets runs at its fixed k on today's
sampler (the move declines at a fixed k), only if it passes A. C6 is per-draw weights on r.

**Length.** Each run is extended until k, the mean count over rows and 90 percent of rows' F reach 100
effective draws (bulk, pooled over chains, batch means of at least 10 batches), after 1,000 sweeps of
burn-in a chain, and is capped at 100,000 kept sweeps. A run that reaches the cap short of its target is
kept, flagged unsettled, and read in the main reading with its flag; a cell with more than 20 percent of its
replications unsettled is reported apart from its band's average; if more than a third of a band's cells
are reported apart, that band's long reading is reported, not read, and the finding says so. By the probes,
narrow and moderate cells take 8,000 to 25,000 kept sweeps, wide cells about 50,000, and the datasets whose
trees mix slowly (solder, mite, and so BCI's like) reach the cap.

**The default reading** is a run of today's sampler at bart()'s defaults and needs no mixing of k: it is
what users get. If the collapsed move later lands for negative binomial, the default reading changes and
the long reading does not; the finding says so.

### Considered: integrating k out over a fixed-k grid

The first version of this plan fitted each replication at a grid of fixed k values and integrated k out by
path sampling, the integrand L - k^2 m SS / A^2 at each node. Answering the critique, it would need: a grid
extended both ways with a k = infinity node (the engine accepts it; the integrand goes to 0 there), since on
narrow and null data most of k's posterior lies above the first version's top of 4; nodes sized for the
single-row estimands, not only the integrand; a log k Jacobian, empty leaves dropped from L and SS, and a
full burn-in at the first node; and honest error. The integrand's Monte Carlo error by batch means is about
1.7 times its spectral estimate, and on roaches, its easiest dataset, two independent chains' first 1,000
draws at k = 1.6 differed by 3.3 combined standard errors. Sized to hold log p(y | k) within 0.5 nat across
the posterior's bulk, it costs about 97 core-hours on the main run's cells against the move's 65, and stays
an approximation; on roaches its k quantiles also disagreed with the move build's ([Probes](#probes)). So
the grid is not used.

## Probes

Run 2026-10-29 on the laptop (arm64 macOS), two cores: about 2.2 core-hours by this plan (1.6 for the first
version, 0.6 for this revision), and about 0.2 by the critique, whose probes are cited as such. Scripts and
outputs are not checked in.

**Real spreads.** S over all rows and over rows with a fitted mean of at least 0.01 from four fits per
dataset ([Real datasets](#real-datasets)): over identifiable rows each dataset's range across fits is
within 7 percent of its midpoint. crimtab is the only dataset with many rows below 0.01 (11 to 33 percent
by fit), and its spread over all rows depends on the fit (2.75 for s(), 3.80 for te()) as its identifiable
spread does not (2.12 to 2.43). A BART fit at k fixed at 1 gave all-row spreads of 1.66 (solder), 1.72
(aids), 1.97 (mite) and 3.42 (crimtab): wider, since it carries interactions and, at a shape floored at 1,
dispersion.

**The prior's spread by candidate.** 2,000 forests drawn from the prior at unit spread (k = 3) on six
designs (crimtab, roaches, solder, aids, mite, and 5 correlated normals at n = 1,000) give S at unit spread,
median 0.48 to 0.67; S under a prior on the spread is that times the spread. They reproduce the first
study's default (median 0.91, 90th percentile 3.1, 99th 14 at crimtab's design). A preview of A1 on these
six designs only; stage 1 gives the full Part A:

| candidate | P(S > 4) |
|---|---|
| C0, chi(1.5, 2) | 0.054 to 0.084 |
| chi(3, 1.25) | 0.009 to 0.020 |
| chi(4, 1.03) | 0.003 to 0.0085 |
| chi(5, 0.9) | 0.001 to 0.004 |
| chi(7, 0.75) | 0.0002 to 0.0014 |
| chi(10, 0.62) | under 0.0004 |
| chi(1.5, 3) | 0.029 to 0.046 |
| C0 with k at least 0.6, 0.75, 1 | 0.0006 to 0.0067; under 0.0024; under 0.0004 |
| half-normal, median 1.57, 1.3 | 0.004 to 0.019; 0.001 to 0.0074 |
| k fixed at 1.5 or 2 | 0 |

The critique's A2 probe, drawing c as log of the mean count with f centered at 0 (not this plan's
convention, so the realized mean runs far above the label) and counting only the top row, found a refused
count in 6 to 9 percent of default datasets at mean counts of 10 and 100, 0.7 to 1.8 percent under chi(4,
1.03) and 0.3 to 1.1 percent under chi(5, 0.9). Under this plan's convention A2 is new for every candidate.

**k under today's sampler.** Seven datasets (roaches, solder, aids, mite, and three simulated with the 5
correlated normals and a binary column at n = 1,000: "wide", S = 2 at a median count of 0.3; "narrow",
S = 0.3 at 5; "zeros", S = 1 at 0.3), 4 chains of 16,000 sweeps under the default: k's autocorrelation
time 440 (narrow), 1,140 to 1,630 (aids, roaches, wide, zeros), 3,720 (solder) and about 8,000 (mite, 7
effective draws). The estimands mixed in 6 to 95 sweeps on four, 590 at wide's single rows, 1,100 to 1,900
on solder and up to 6,100 on mite; a default-length run's 90 percent intervals held 0.85 to 0.92 of the
long run's draws on the first four, 0.69 to 0.92 on wide, 0.73 to 0.79 on solder, 0.54 to 0.92 on mite.
chi(5, 0.9) did not reliably shorten k's autocorrelation. The critique's runs on Traffic, Insurance and null
data at n = 500 put the default's median k at 4.7, 5.7 and 6.1, with 73, 97 and 98 percent above 4.

**The scratch move build.** Built from collapsed-move-main 0a9e9e41 with the negative binomial evaluator and
predicate (about 15 lines), installed in its own library. The same runs under it, and Traffic (4 chains of
1,000 burn-in and 20,000 kept) under both builds:

| data | k's autocorrelation, today / move | slowest estimand's | k's 5th, 50th, 95th: today | move | ms a sweep, today / move |
|---|---|---|---|---|---|
| Traffic | 620 / 37 | - | 2.97, 4.90, 7.12 | 3.08, 4.95, 7.10 | - |
| narrow | 441 / 51 | 36 / 8 | 2.29, 2.98, 3.81 | 2.30, 2.98, 3.85 | 2.8 / 3.3 |
| roaches | 1,200 / 236 | 39 / 25 | 0.92, 1.31, 1.82 | 0.82, 1.22, 1.72 | 1.5 / 1.7 |
| wide | 1,550 / 506 | 590 / 513 | 0.61, 0.72, 0.86 | 0.58, 0.70, 0.83 | 5.2 / 6.1 |
| solder | 3,720 / 1,430 | 1,910 / 2,390 | 1.05, 1.28, 1.56 | 0.82, 1.10, 1.44 | 8.1 / 10.0 |
| mite | 8,000 / 1,600 | 6,110 / 4,020 | 0.52, 0.64, 0.79 | 0.54, 0.63, 0.75 | 3.6 / 4.3 |

Where today's sampler converged (Traffic with 129 effective draws of k, narrow with 136) the two agree within
Monte Carlo error, in k and in the estimands. On roaches and solder they differ in k's lower tail, where
today's chains had 50 and 17 effective draws and the move's 255 and 42; on solder the move sides with the
grid of the first version (median 1.13). Corrected 2026-10-29: the critique's long runs on roaches settle
it for the move build ([Amendments before results](#amendments-before-results)); the disagreement came from
today's short runs and the first version's grid. The estimands' intervals agreed on every dataset where
both runs converged. The move costs 17 to 24 percent more a sweep, and on solder, mite and wide's single
rows the trees, not k, set the pace.

**The grid of the first version**, on five datasets (two path chains, step 0.1 in log k, 100 sweeps of
burn-in after each change, 1,000 kept a chain per node): k's 5th, 50th and 95th percentiles 0.95, 1.26, 1.83
on roaches (move build 0.82, 1.22, 1.72), 1.43, 1.75, 2.09 on aids, 0.60, 0.71, 0.84 on wide, 0.96, 1.13,
1.38 on solder and 0.57, 0.68, 0.78 on mite; the integrand's effective draws a node were 4 to 220 on roaches
and aids and 3 to 50 on the rest. The units check passed: under the prior at fixed k the integrand's mean was
-0.50, -0.39 and -0.88, standard errors about 0.5, at k 0.8, 2 and 4.

**Cost.** Milliseconds a sweep of today's sampler, fitted to the seven direct runs, are about 0.16
microseconds times the sum over rows of y + r, plus 0.9 microseconds a row, at 75 trees. r is the drawn
shape, so Poisson-like data cost through n r: solder's 8 ms is 900 rows at a shape near 46.

**A package quirk met on the way**: `dbarts()` given `offset.test = v` with `v` an existing variable whose
value is NULL errors "cannot find test offset 'v'"; the probes omitted the argument instead. Reported, not
part of this study.

## Cost

Projected from the probes' cost model and their autocorrelation times, the move build at 1.2 times today's
cost a sweep; the box's times are assumed equal to the laptop's. The pre-registration replaces these with
stage 2's full-length runs, which cover every band.

| stage | where | core-hours |
|---|---|---|
| 0 harness, truths, checks (the move build's SBC of k included) | laptop | about 2 |
| 1 prior predictive, 20,000 datasets a cell | laptop | about 1 |
| 2 the long reading's check, full length, every band | laptop | about 4 |
| 3 smoke | box | about 1 |
| 6 main run: long runs 65, real-outcome folds 8, C0's default readings 2 | laptop about 23 (narrow, middle, small-sample wide); box about 52 (wide, folds) | about 76 |
| 6, if C1 passes A: its fixed-k runs | box | 28 more |
| 7 guard: one candidate's default readings (a second only if the first fails; C4 or C5 on a scratch build) | box | about 2 (4) |
| 8 confirmation, only if called: a third of the main run, with its guard | box | about 25 (35) |
| expected: C1 failing A, no confirmation, one guard | | about 86 |
| at most | | about 150 |

The main run: about 966 replications. Narrow cells 20 replications each, about 6 core-hours; middle 10,
about 5; small-sample wide 20, about 12; wide 10 (BCI 4), about 43 by the long runs, the slow-tree
datasets at the cap the largest share. Candidates whose weights fail need their own runs: assumed for four
candidates on every narrow and small-sample replication and for one in two elsewhere.

Priced against the alternatives on the same cells: today's sampler's long runs, to 100 effective draws of k
at its autocorrelation times, about 181 core-hours; today's sampler where k mixes fast and the grid
elsewhere, sized honestly, about 97 and approximate; the move build, about 65. The move build is chosen.

**Stop**: 160 core-hours over all stages, counted from the harness's timing logs; at the stop the reading is
taken on what has run and the finding says so. If stage 2's projection exceeds 150, cuts are made in this
order before the pre-registration commit, sparing the narrow and small-sample cells, where candidates can
differ: the BART-fitted truths' wide-band cells to 5 replications; BCI's cells; the slow-tree datasets' cap
from 100,000 to 50,000 kept sweeps; the wide band's replications from 10 to 6.

**Where and when.** Stages 0 to 2, about 7 core-hours, run on the laptop's two cores, about 4 hours of wall
time, inside the roughly 20 hours the box needs for the k study's first tranche; stage 2 waits on the
maintainer's ruling after stage 1. The smoke and the main run's wide cells and folds, about 53 core-hours,
go to the box's four cores after that tranche, about 13 hours; the laptop takes the narrow, middle and
small-sample cells, about 23, about 12 hours on two cores. A confirmation adds about 7 hours on the box.

## Stages and gates

Each gate is a stop: a failure is a harness defect to fix, not a finding, unless it says otherwise.

0. **Harness and checks** (laptop). Data export, truths, the harness, the scratch builds, and these checks:
   - **The move build samples the posterior**: SBC of k on 200 datasets drawn from the default prior at
     n = 50 (ranks uniform by a chi-squared test at 0.01), and on two of them k's draws match today's
     sampler's within Monte Carlo error.
   - **Reweighting is exact**: on three prior-predictive datasets, chi(5, 0.9)'s posterior reweighted from
     C0's run matches its own run's (k's quantiles and F's intervals within two combined standard errors).
   - **The prior draw is the fitted model's prior**: cut points, center c and s = A / k as bart() sets them
     (the first study's checks, reused); Part A's convention reproduces a fit's c; the top-rows shortcut for
     the largest count agrees with full draws on 1,000 datasets a cell.
   - **Truth-blind**: no run length, weight, interval, group definition or acceptance reads the truth, and
     the harness's scoring against the truth is off until the pre-registration commit turns it on.
1. **The prior predictive** (laptop): every candidate on every A cell and real design; Part A's full table
   (A1, A2 and A3 for every candidate, with S's and S_id's quantiles) written to the harness output. No
   fit.
   **The maintainer rules on the criterion** with that table in hand (open call 1, and the others), before
   any Part B fit. The ruling is recorded as a ledger entry quoting the maintainer.
2. **The long reading's check** (laptop): 16 replications across every band (the real datasets but BCI,
   Insurance and Traffic, as fitted; three widened; D4, D5, D8 and two small-sample cells), each at full
   length on the move build, and on the narrow and middle ones also a long run of today's sampler to 100
   effective draws of k. Passes when, where both converge, k's 5th, 50th and 95th percentiles and F's
   interval endpoints agree within two combined Monte Carlo standard errors (both runs' errors by batch
   means). It also measures each band's run lengths, weight failures, unsettled share and times, from which
   the pre-registration's counts and costs are set. Blind: the harness reports agreement, lengths and times
   only, candidates under coded labels. If the move build fails and cannot be fixed (it is a harness
   defect until shown otherwise), the study stops and goes back to the maintainer with a costed choice
   between today's sampler's long runs on fewer cells and the grid.
3. **Smoke** (box): every cell, one replication, a tenth of every length, every run kind and the
   real-outcome folds. Blind, as stage 2; it checks the pipeline end to end, not lengths or costs.
4. **Review** (opus) of the harness, the scratch builds and the smoke; fixes.
5. **Pre-registration commit** (below), before any study fit is scored.
6. **Main run** (box and laptop, as [Cost](#cost) splits it): every cell's replications, their long runs
   and C0's default reading; the real-outcome folds. At 10 replications a cell the reading is taken once:
   if every candidate's step 2 verdict clears its bound by at least the bar's half and step 3's leader
   leads the next untied candidate by at least 0.01, the cell's remaining replications are skipped (the
   house's rule to stop once a result is settled).
7. **Reading and guard**: steps 1 to 4 of [The decision rule](#the-decision-rule); then the recommended
   candidate's default reading on the main run's replications (C4 or C5 on a scratch build carrying only its
   k draw), and the next one's only if it fails; step 5.
8. **Confirmation**, only if step 6 of the rule calls it.
9. **Finding**, its opus review, a ledger entry for the calls made running it, and the result to the
   maintainer ([How the result goes to the maintainer](#how-the-result-goes-to-the-maintainer)).

## Pre-registration

Made 2026-10-29 by the commit that turns the harness's scoring on, with this plan, the harness, the two
patches and the box script at their state in it; no study fit had been scored.

Committed to this plan, with the harness at that commit, after stages 0 to 4 and before any study fit is
scored: the maintainer's ruling on the criterion, quoted; the decision rule; the candidates and their exact
priors, C6's weights and C7's formula; the run-length rule, cap and unsettled rule; the replication counts
and the truths' formulas and seeds; stage 1's Part A table and verdicts; stage 2's result and the per-band
costs and projected total from it. Fixed by this plan and not changed by stages 0 to 4: the bands, bars and
bounds, the candidate list, the cells, the stop and the order of cuts. What stages 0 to 4 may change is
listed, and only that: replication counts (down, to fit the stop, by the order of cuts), the weight-failure
threshold (to keep own runs within cost), and a widened truth's center (by the stated rule). Any other
change after a fit is scored is a deviation, disclosed in the finding.

Disclosed with the pre-registration, from the review of the Part B fixes:

- **The cell rule reads high near the bar.** A cell's bound linearizes its gap difference and leaves out
  C0's folded-gap noise near nominal, so the rule fires at the bar, not below it: in the reviewer's
  constructed case on F, a true difference of 0.021 flagged 14 of 68 cells, and 0.010 and 0.015 none.
- **The information fractions assume every cell gets both blocks**, so a look's fraction is understated
  when fewer run (0.354 at look 1 in the reviewer's case), which raises the early boundaries; the error
  spent still sums to 0.05, the final look's fraction being 1. Each look recomputes the earlier looks'
  boundaries from its updated variances. The cell rule and the width guard are fixed-level t tests
  repeated at each look, outside the spending.
- **Prices not measured**: C4 at its floor and C5 on null data have not run at full length; their own runs
  are priced at C0's length.
- **Paths not coded**: step 4 and step 3's first tie-break, which stage 1's results rule out (the default
  fails A, and no candidate passes it with margin). The harness records both results as constants, and
  the decision stops with an error if either is changed.

## Amendments before results

Made 2026-10-29 on the orchestrator's calls after the critique's re-check of this plan (verdict: build
stages 0 and 1), before stage 0 ran and before any fit was scored. Where they differ from the text above,
they govern.

- **The review before stage 0.** The critic's reading of the move's negative binomial lines and its long
  runs on roaches count as the review of the scratch move build before stage 0; the stage 4 opus review
  stays. Its runs, 2 chains of 200,000 sweeps on today's sampler and 2 of 60,000 on the move build, put k's
  5th, 50th and 95th percentiles in agreement within batch-means error, so the move build's lower tail on
  roaches is the right one. The disagreement [Probes](#probes) reported there came from today's short runs
  (50 effective draws of k) and the first version's grid.
- **Reweighting, a second check.** Beside the 30 percent rule, Pareto-smoothed importance sampling: a
  candidate's weights on a replication are used only when the fitted Pareto shape of their upper tail,
  k-hat, is at most 0.7 and they keep 30 percent of the run's effective draws; otherwise that candidate gets
  its own run on that replication. Candidates whose weight against C0, pi_c(k) / pi_0(k), grows without
  bound as k grows get their own runs, unweighted, wherever the likelihood is flat at large k (the narrow
  band, the null cell, and any replication where C0's run puts more than 5 percent of its draws above
  k = 4): chi(1.5, 3); both half-normals, whose density in k falls as 1 / k^2 against C0's exp(-k^2 / 8);
  and C7 wherever its chi scale exceeds C0's 2, which on narrow data it does.
- **Own runs for C4 and C5, stated precisely.** They use a scratch build that changes two things against
  the move build. k's own step: for C4, k^2 drawn from its conjugate gamma truncated to k at least the
  floor, by inversion; for C5, a slice step on log k. And the move's k term: the move's log density in
  v = log r is log pi(k e^-v) - v + sum_i l_i(e^v f_i + o_i), with pi the density of k, of which the
  shipped -nu v - (k e^-v / s)^2 / 2 is the chi(nu, s) case up to a constant. The scratch build puts the
  candidate's own pi there: for C4, chi(1.5, 2)'s, and minus infinity where k e^-v falls below the floor
  (the slice, started at v = 0 inside the support, never leaves it); for C5, log 3 - 2 log k -
  (3 / k)^2 / (2 sigma^2), sigma the half-normal's scale (1.93 for median 1.3, 2.33 for 1.57). Each such
  build is made only if its candidate passes A, and passes stage 0's drift test before it is used.
- **Stage 0 gains two checks.**
  - The reweighting check adds chi(1.5, 3) and C7 on null data (D5's design): their k-hat against C0's
    run, and their own runs against the reweighted reading where k-hat is at most 0.7.
  - A drift test from the truth, with offsets, as the k study ran one: per replication, 80 rows of 3
    correlated normal columns and an exposure offset log-uniform over 1 to 100; trees, k0 from chi(1.5, 2),
    leaves at k0 and a shape r0 from the grid prior, held fixed in the sampler so that the start is exact
    (the shape's own draw is shipped code; the move's refresh then redraws the latents only); counts from
    the model at the build's center; the chain started at that state, an exact posterior draw once the
    latents are drawn given the counts. log k averaged over lags 100 to 1,000 within a replication must
    not drift from log k0: |z| at most 3 over 2,000 replications. It must also have power: a mutant of the
    move whose Jacobian carries one power of r too many (-(nu + 1) v for -nu v) must give |z| above 5, or
    the test does not count as passed.
- **Budget.** About 100 to 150 core-hours expected, stop at 160. The candidates likely to pass A are the
  thin ones (chi with 4 or more degrees of freedom, C4, the narrower half-normal), and those need their
  own runs on wide cells more often than one replication in two: C4 wherever k's posterior sits at or
  below its floor, the half-normal and C7 wherever the likelihood is flat at large k. The cap of 100,000
  kept sweeps is per run, own runs included.
- **The width guard** (10 percent on each narrow-band cell) may exclude the thinner chi priors: on null
  data the critique's probe found their intervals about 20 percent wider than the default's. The summary
  says so.

## Stages 0 and 1

Run 2026-10-29 on the laptop's two cores by the harness's stage 0 and 1 actions, on the move build and on
main 8a78507c ("today"), about 1.5 core-hours in all, against the 3 budgeted. No fit was scored against a
truth.

**Sized** on the maintainer's instruction, the same day, that runs and replications be sized to the
question each answers and no further:

- The drift test ran 1,000 replications on each build, in batches, not 2,000. The statistic's spread across
  replications is 0.29 and the mutant's error moves log k by only about 0.05 within the window, so |z| above
  5 on the mutant takes about 1,000; the move build ran to the same count, at which a drift of 0.04 in log
  k, under a third of k's narrowest posterior sd on the real data, shows at |z| 3 with 90 percent power.
- The checks' prior datasets are conditioned on their total count: 50,000 in the drift test, above every
  fitted real dataset's, and 10,000 in the SBC. A condition on the counts alone keeps each check exact. It
  bounds the cost, which grows with the total: without it, single replications ran for many minutes.
- The SBC keeps the plan's 200 datasets, at which a bias of a quarter of k's posterior sd moves the mean
  rank by 3.4 standard errors. Each chain runs until 99 draws spaced at twice k's autocorrelation time
  fit, capped at 100,000 sweeps, not a fixed length.
- Stage 1 drew forests in blocks of 2,500 a design and stopped once every verdict was resolved (the
  estimate two standard errors from its bound, or its standard error under 5 percent of the bound): 5,000
  forests on 23 designs, 20,000 on mite and on crimtab (rerun there for the options below).
- The BART-fitted truths are 2 chains of 1,000 kept draws, not long fits: a truth is a fixed function once
  stored, so its Monte Carlo error does not enter the study.

**Stage 0.**

- *Drift test from the truth.* Move build, 1,000 replications: drift of log k over lags 100 to 1,000
  +0.008 (se 0.009), z +0.89; the shape held at r0 in all. The mutant: +0.050 (se 0.009), z +5.47. Passes,
  with power.
- *SBC of k.* 200 datasets at n = 50: ranks' chi-squared 7.1 on 9, p 0.63; mean rank z -0.92. Passes.
  k's autocorrelation time had median 28 sweeps and 90th percentile 171; 9 datasets, with totals in the
  thousands, reached the cap at times of 520 to 2,300 sweeps, and without them p is 0.80. The shape's
  ranks, reported and not gated, give p 0.015, from one decile (37 of 200) filled by datasets drawn at
  shapes 12 to 20. The shape's draw is shipped code the move does not touch, and the chains were thinned
  by k's mixing, not the shape's; stage 4's review checks it. 0.5 core-hours.
- *The builds on two datasets.* SBC datasets 1 and 2, 2 chains of 100,000 sweeps on each build: k's 5th,
  50th and 95th percentiles within 0.63 combined standard errors (move 2.05, 3.61, 5.80 against today's
  2.01, 3.63, 5.81; and 1.10, 1.58, 2.21 against 1.11, 1.59, 2.22). Passes. The move build's effective
  draws of k were 13,700 and 2,400 against today's 1,100 and 540.
- *Reweighting.* chi(5, 0.9) from C0's run on three prior datasets (200 rows): weights kept 37 to 95 percent,
  k-hat at most -0.81; k's 5th, 50th and 95th percentiles reweighted and from its own run within 2.1
  combined standard errors (the nine z's give chi-squared 14.6 on 9, p 0.10), rows' interval endpoints
  within two in 92, 94 and 97 percent of rows. Passes. On null data (D5's design, two replications)
  chi(1.5, 3) kept 11 and 20 percent with k-hat 0.48 and 0.54, C7 1 and 4 percent with k-hat 0.79 and 0.85;
  the reweighted 95th percentiles of k differed from their own runs' by 6 to 24 standard errors. So
  reweighting fails there, as amended, and the 30 percent rule, not k-hat, is the one that catches
  chi(1.5, 3). Own runs stand as amended.
- *The prior draw.* A = 3, and at k = 3 the rows' prior variance is 0.93 to 1.02 on four designs (1
  expected). Cut points are identical to bart()'s, and the sampler's c is the formula's. A fit's c on Part
  A's datasets is centered on the log of the cell's mean count (medians within 0.002), its central 90
  percent within 0.11 at mean counts of 1 and above and within 0.3 at 0.1. The top-rows shortcut for a
  count above the cap is within 1.3 standard errors of full count draws, including the one checked cell
  where a count passes the cap (10,000 rows, mean 100: 0.67 against 0.80 percent).
- *Data and truths.* 13 datasets exported with their additive truths; the four fits' spreads reproduce
  [Real datasets](#real-datasets). The BART truths at k = 1 have S_id 0.29 (Insurance) to 2.33
  (crimtab).
- *What stage 0 changed in the harness.* Effective draws are the smaller of batch means (10 batches) and
  Geyer's estimate: batch means alone read high early in a run, where the rule checks first. The
  simulated offset design's exposure is scaled to mean 1, so a cell's label is a row's mean count as in
  the other designs. Unscaled, a mean of 100 put rows at 10^5 by exposure alone.

**Stage 1, Part A.** Percent of prior datasets: A1 and A2 at the worst design (at most 1 to pass), A3 at
crimtab, the lowest real dataset for every prior but the two below the line (at least 1 to pass), read at
crimtab's widest fit (2.43, the rule as written) and at the middle of its four fits' range (2.27).

| prior | A1 | A2 | A3 at 2.43 | A3 at 2.27 |
|---|---|---|---|---|
| C0, chi(1.5, 2) | 8.41 | 0.77 | 5.52 | 8.90 |
| C2 chi(3, 1.25) | 2.00 | 0.06 | 1.46 | 3.18 |
| C2 chi(4, 1.03) | 0.99 | 0.04 | 0.65 | 1.71 |
| C2 chi(5, 0.9) | 0.52 | 0 | 0.39 | 1.15 |
| C2 chi(7, 0.75) | 0.17 | 0 | 0.12 | 0.57 |
| C2 chi(10, 0.62) | 0.04 | 0 | 0.07 | 0.26 |
| C3 chi(1.5, 3) | 4.62 | 0.43 | 3.43 | 5.02 |
| C4, k at least 0.6 | 0.89 | 0 | 0.68 | 2.14 |
| C4, k at least 0.75 | 0.27 | 0 | 0.33 | 1.15 |
| C4, k at least 1 | 0.02 | 0 | 0.07 | 0.33 |
| C5, median 1.3 | 0.95 | 0 | 0.80 | 2.07 |
| C5, median 1.57 | 2.29 | 0.01 | 1.54 | 3.96 |
| C1, k = 1.5 | 0 | 0 | 0.04 | 0.19 |
| C1, k = 2 (also below 1 on aids, BCI, mite) | 0 | 0 | 0 | 0 |
| C7, at each real design's own spread | 7.72 | 0 | 0.01 | 0.02 |

Standard errors near a bound are about 0.03 points on A1 and under 0.1 on A3. C6's shape weights move A2 by
at most 0.2 points. At the simulated spreads C7's A1 is 0, 0.69 and 8.41 percent (0.25, 1 and 2).

- The default fails A1 at every design, 4.0 to 8.4 percent; its A1 is largest at the real designs.
- A2 separates nothing: under the convention every prior refuses a count in under 1 percent of datasets,
  the default's worst 0.77 percent (10,000 rows, mean count 100). The critique's 6 to 9 percent came from
  a center not set by the data.
- A1 and A3 pull apart at crimtab, and as written no candidate passes A. Every candidate that passes A1
  leaves crimtab beyond its 99th percentile, and every one that keeps crimtab inside it fails A1. The
  pinch is crimtab's widest fit: the A1-passing candidates' 99th percentile of S_id there is at most
  2.39, against 2.43. Read at the middle of crimtab's range, five pass: chi(4, 1.03), chi(5, 0.9), C4 at
  0.6 and 0.75, and C5 at 1.3. None passes with margin. Every other real dataset lies inside the 99th
  percentile of every A1-passing prior but k fixed at 2.
- C7 fails both ways: scaled to a main-effects Poisson fit, its spread leaves crimtab beyond its 99th
  percentile and exceeds A1's bound on five real designs (BCI the worst, 7.7 percent).
- Migration, reported only, lies beyond the 99th percentile only under k fixed at 2.

**Proposed sizes for Part B**, for the maintainer's ruling with the criterion (open call 8), each sized to
the smallest difference that would change a verdict: the tie margin, 0.005 in a class's coverage gap, and
half the bar, 0.01, for a confirmation. Until stage 2 measures it, a paired coverage difference's spread
across replications is assumed 0.03 a cell in the middle and wide bands and 0.10 in the narrow.

- Replications: narrow and small-sample cells 20, as planned; middle and wide cells 5 (BCI 2), not 10. A
  class's gap then has a standard error near 0.0024 against 0.002, still under half the tie margin, and the
  wide band, where candidates differ least, costs half. Saves about 24 core-hours.
- Effective draws: 100 for k and the mean, which the weights need; 50 for the rows' interval endpoints,
  whose Monte Carlo error is shared by every candidate reweighted from one run and so mostly cancels in
  their differences. This halves the slow-tree datasets' run lengths and saves about 12.
- Own runs: only for candidates passing A and only where the weights fail, which stage 0 found on null data
  and not on prior data; about one in ten middle and wide replications, not one in two. Saves about 13.
- Confirmation: on close bands only, in blocks of 5 replications until the margin's interval clears 0.01,
  in place of a third of the main run on every cell; about 6 core-hours against 25.
- Kept: the cells, the cap of 100,000 sweeps (mite needs it), and the interim stop at 10 replications.

With five candidates passing A (open call 1's recommendation), the study comes to about 50 core-hours at
these sizes against 85 at the plan's; stage 2 re-prices both.

## Amendment after stage 1

Made 2026-10-29 after stage 1, on the maintainer's ruling of the criterion and on the orchestrator's calls
after the critic's audit of the proposed Part B sizes against the maintainer's ask that runs be sized to
their questions; before stage 2 ran and before any fit was scored. Where they differ from the text above,
the proposed sizes for Part B in [Stages 0 and 1](#stages-0-and-1) included, they govern.

- **The criterion, dec-B502.** The maintainer: "OK, use your recommendation for this item." A3 reads each
  real dataset's S_id at the middle of its range across the four fits, at 1 in 100; A1 and A2 as written;
  Part B as recommended. Five candidates pass A and go to the fits: chi(4, 1.03), chi(5, 0.9), C4 with k at
  least 0.6 and at least 0.75, and C5 at median 1.3. None passes with margin, so a tie goes to the smaller
  code change. crimtab is among the fitted cells. The candidates that fail A are reported on B where their
  weights hold and cost nothing more; C1 is not run.
- **Replications.** Narrow and small-sample cells 20; middle 5; wide 5, except crimtab's cells, as fitted
  and widened, 10. The interim stop at 10 of 20 replications is dropped: under O'Brien-Fleming boundaries
  it needs a lead of about 0.015 and would rarely fire.
- **Effective draws, per candidate.** The target of 100 effective draws for k, the mean and 90 percent of
  rows is counted per candidate after reweighting: a run is extended until every candidate whose weights
  hold reaches it, or the cap. A candidate below 50 after reweighting on a replication gets its own run
  there, beside the Kish and k-hat rules. Each candidate is scored against its own expected coverage,
  0.9 N_c / (N_c + 1), N_c its effective draws on the row, in place of 0.9: scored against 0.9, a
  reweighted candidate at 15 to 25 effective draws would sit 0.035 to 0.056 short, a spurious loss of up
  to twice the bar. Stage 2 checks the correction on its long runs: an interval from two of a run's
  chains, reweighted, covers the other two chains' draws at that rate.
- **A guard on the widest cells**, pre-registered: in each class, the average gap over crimtab's cells and
  the cells widened to S_id 3 (about 13 cells) is not worse than the default's by the bar at its upper
  bound. The tail cuts and the thin priors put their risk there, and the band average dilutes it.
- **Confirmation, group-sequential.** The main run is the first look, followed by at most two blocks of 5
  replications. Boundaries are O'Brien-Fleming's, z of 3.47, 2.45 and 2.00 at the three looks, by
  Lan-DeMets spending at the realized information fractions if the blocks are unequal. The thresholds
  tested are the decision targets: 0.02 for non-inferiority, 0.005 for the tie. A result unresolved at the
  last look reads as a tie for step 3 and a fail for step 2. This replaces step 6's third of the main run.
- **Cells and folds.** BCI's fitted cells are dropped (about 5 core-hours); the real-outcome folds run one
  repetition (about 4).
- **Budget.** About 62 core-hours expected, re-priced from stage 2. The stop comes down from 160 to about
  110 unless stage 2's prices say otherwise.
- **Revisit after stage 2.** If the measured paired spread on narrow and small-sample cells exceeds about
  0.12, those cells need more replications, and stage 2's record says so. Stage 2 measures it blind, as
  blinded sample-size re-estimation does: on 5 replications each of three narrow and four small-sample
  cells, the paired difference in F's 90 percent coverage, candidate less default, is computed against
  the truth under coded labels, and only its spread within cells, centered and pooled, is summarized; no
  mean or direction is printed. Candidates whose weights fail there are left out of it where their own
  run needs a scratch build.

## Stage 2

Run 2026-10-29 on the laptop's two cores by the harness's stage2 action, under the
[Amendment after stage 1](#amendment-after-stage-1): about 2.2 core-hours, against the 4 budgeted. The
move build ran C0's long runs on 18 cells at one replication (the stage's list: nine real datasets as
fitted, roaches and crimtab widened to S_id 3 and quakes to 2, D4 at 3, D5, D8, and two small-sample cells;
and aids at 30 rows for the spread), the seven spread cells to five replications, and the chi candidates'
own runs where their weights failed; today's sampler ran C0's long runs on the narrow and middle cells but
solder. No fit was scored against a truth but the blind spread.

- **The long reading's check passes.** Where both builds settled (seven cells, all of them), k's 5th, 50th
  and 95th percentiles agree: 21 z's, chi-squared 28.0 on 21, p 0.14, the largest |z| 2.55 (grouseticks'
  95th); the rows' interval endpoints within two combined standard errors in 92 to 97 percent of rows.
- **Lengths.** All 46 runs settled. C0's own target took 2,500 to 5,000 sweeps a chain on the narrow and
  middle cells (solder 5,000), 2,500 to 15,000 on the wide (mite 7,500; neither slow-tree dataset neared
  the cap the probes projected), and up to 20,000 on the small-sample. Counting the target per candidate
  after reweighting (the amendment) lengthened the widened-to-3 and small-sample runs 1.7 to 4.7 times
  (crimtab widened to 3: 7,500 to 22,500), where a candidate's weights held with a Kish share near 0.3 to
  0.5.
- **Weights.** Own runs a replication: narrow 2.4 (chi(4, 1.03) on 67 percent, chi(5, 0.9) on 73, the
  half-normal on all, its weight unbounded where the likelihood is flat); middle none; wide 0.2 (C4 at 0.75
  on crimtab widened to 3, whose k sits at 0.45 to 0.70); small-sample 1.15 (C4 at 0.75 on 40 percent, C4 at
  0.6, chi(5, 0.9) and the half-normal on 25 each). C4's and C5's own runs need their scratch builds.
- **The expected-coverage correction.** An interval from two of a run's chains covered the other two's
  draws at 0.9 N / (N + 1) within about 0.001 at 100 or more effective draws; at 50 to 100 the formula
  under-predicts the coverage by 0.004 to 0.006 (0.894 against 0.888, reweighted candidates over 4,610
  rows), and by 0.02 to 0.03 at 25 to 50: N as counted, the smaller of batch means and Geyer's estimate,
  is conservative for the interval's tails. Before the pre-registration N_c is either read from a tail
  estimate of effective draws or the correction calibrated on these runs.
- **Paired spread, blind** (a candidate's F coverage less C0's, sd across replications within a cell,
  pooled; candidates under coded labels, no mean read): narrow, three cells, 0.004 on 48 df, the largest
  candidate's 0.009; small-sample, four cells, 0.019 on 65 df, the largest 0.045. Both are far under the
  0.12 that would call for more replications; at 0.009 the narrow band's 7 cells resolve the tie margin
  at 5 replications.
- **The small-sample cells as written are degenerate.** Widened to S_id 3 at a mean count of 5, only 8 to
  30 percent of roaches', epil's, aids' and quakes' rows keep a mean of 0.01, and a subsample of 30 is
  mostly zeros (a stage 2 smoke drew a dataset of all zeros). Stage 2 widened them by the widened cells'
  rule instead (the dataset's center, lowered for the caps) at shape 3. D4 at S 3 keeps 17 percent of its
  rows; its k sits at 0.06 to 0.5, and its runs were the longest relative to their size.
- **Re-priced**, from these runs (the move build's seconds a sweep fitted on rows and counts; a band's
  sweeps a replication, own runs included, C4's and C5's priced at C0's length), at the amended sizes
  (760 replications): long runs about 50 core-hours (narrow 10, middle 2, wide 15, small-sample 23), C0's
  default readings about 3, the folds 4, the guard 2, the smoke 1: about 60 expected before confirmation.
  A block of 5 replications costs about 2 on the narrow or middle band, 6 on the small-sample and 14 on
  the wide; both blocks on every band come to about 110, the amended stop.
  Two cuts the measurements support: replications on the narrow band at 5, not 20 (saves about 7), and no
  extension for the candidates, a candidate below 50 effective draws at C0's stop getting its own run
  (saves about 7).
- **Before the main run**: the Part B scoring and the folds in the harness; C4's and C5's scratch builds,
  each passing stage 0's drift test; the small-sample cells' definition, the correction's N and the two
  cuts above settled; then the smoke, the review and the pre-registration commit.

## Amendment after stage 2

Made 2026-10-29 on the orchestrator's rulings on stage 2's report, before any fit was scored. Where they
differ from the text above, they govern.

- **Small-sample cells.** Each is widened by the widened cells' own rule at shape 3: the dataset's center,
  lowered only as far as the caps on the expected total and the largest mean require, in place of a mean
  count of 5. As written, only 8 to 30 percent of rows kept a mean of at least 0.01 and one stage 2 dataset
  was all zeros: the definition was degenerate.
- **Narrow cells at 5 replications**, not 20 (the first cut). The narrow band's paired spread, 0.009,
  resolves the tie margin at 5.
- **No extension for the candidates** (the second cut). C0's run stops at its own target. A fitted
  candidate whose weights fail, or which keeps fewer than 100 effective draws after reweighting, gets its
  own run. Effective draws are counted as the run rule counts them (the smaller of batch means and Geyer's
  estimate, times the weights' Kish share), in the own-run rule and in each candidate's expected coverage
  alike.

  The ruling took the cut only if the expected-coverage correction is accurate to 0.002 at the own-run
  threshold, by a tail count of effective draws (the endpoints' indicators) or by a threshold of 100,
  whichever costs less. Measured on 23 runs of 17 cells without extension, on 200 rows a run: an interval
  from chains 1 and 2 (all, half and a quarter of their draws) against the other two chains' draws.
  Coverage less 0.9 N / (N + 1), with its standard error over runs:

  | effective draws | bulk count, C0 | bulk count, reweighted | tail count, C0 | tail count, reweighted |
  |---|---|---|---|---|
  | 50 to 75 | +0.005 (0.004) | +0.005 (0.004) | one row | +0.001 (0.008) |
  | 75 to 100 | -0.002 (0.003) | -0.000 (0.002) | +0.014 (0.005) | +0.003 (0.006) |
  | 100 to 150 | -0.001 (0.002) | -0.001 (0.002) | +0.004 (0.004) | +0.001 (0.002) |
  | 150 to 200 | -0.001 (0.001) | -0.001 (0.002) | -0.002 (0.002) | -0.002 (0.002) |
  | above 200 | -0.0015 (0.0009) | -0.0018 (0.0008) | -0.0031 (0.0011) | -0.0035 (0.0009) |

  The tail count does not hold the error under 0.002: it overstates the effective draws of the rows that
  have the most, where most rows are, and on C0's own reading it understates them at 75 to 100. The bulk
  count holds it from 75 up, within 0.0018 in every bin. So the threshold of 100 on the bulk count is
  taken, although it gives more own runs a replication than the tail count at 50 would (narrow 2.25
  against 2.0, middle 0.2 against 0, wide 1.2 against 0.5, small-sample 1.9 against 1.4, of the five
  fitted candidates).
- **Where the main run goes.** The box, after the k study's first tranche, about 10 to 12 hours of wall
  time from 10:30 box time on 2026-10-29. No main run on the laptop.

## Box readiness

Run 2026-10-29 on the laptop's two cores, before any fit was scored, and revised the same day on the Part B
review's findings (its verdict: fix first) and the orchestrator's calls on them.

**C4's and C5's scratch build.** One build serves both candidates and the guard: the move build, with
its k step and the move's k term taking the prior an environment setting names, read at each use and unset
for the shipped chi(1.5, 2). With the move switched off it is today's sampler carrying only the k step,
the build the guard's default readings need (stage 7). Under C4, k^2 is drawn exactly from its gamma
truncated at the floor: by rejection while the floor keeps at least 5 percent of the conditional's mass,
otherwise from an exponential proposal at the floor, the target's ratio to which is largest there. That
departs from the inversion the [Amendments before results](#amendments-before-results) state: both are
exact at a constant expected cost, and inversion stays accurate even far in the tail (checked to a log
tail probability of -839), so the departure is one of code, not of the draw's distribution. Under C5, a
slice step on log k. The move's k term is the candidate's own log pi(k e^-v) - v, as that amendment states.

- *Identity.* With its prior unset, the build draws as the move build, draw for draw; with the move off,
  as today's sampler; and under a floor that cannot bind, as either. One fit: 200 rows, 2 chains of 600
  sweeps.
- *The truncated draw*, ported to R line for line, matches exact inversion on both branches: a
  Kolmogorov-Smirnov p of 0.54 to 0.94 on six cases, the tail branch down to a bound at 14 times the
  conditional's mean.
- *Drift test from the truth*, 1,000 replications each, k0 from the candidate's own prior: C4 at 0.6,
  z -0.11; C4 at 0.75, z +0.80; C5 at median 1.3, z +0.33, its standard error twice the others' since k0
  spreads wider under the half-normal. A stress run at a floor of 1.5 gives z -0.59; there the chain came
  within 10 percent of the floor in 33 percent of replications, against 6 and 10 percent at 0.6 and 0.75.
  Each passes.
- *What the drift tests cannot see.* They have no power on C4's tail branch: the review's mutant, which
  accepts every exponential proposal there, gave z -0.05 at the stress floor and went undetected, since
  the chain takes that branch too rarely to move log k. That branch's exactness rests on the analytic
  argument (the target's ratio to the proposal is largest at the floor, so accepting with that ratio
  over its value there is exact) and on the R port's match to inversion on that branch.
- *The base.* The two patches are checked in beside the harness, each stating its base: the move patch
  against main at 8a78507c, which pins "today" (main has since moved), and the kprior patch on top of it.
  Applied to a fresh tree of that commit they give the laptop's libraries' source exactly.

**On the box, before any job**, by the checked-in box script: copy the exported data with its truths
(the counterfactual truths included); build the three libraries (today, move and kprior) from a fresh
tree of 8a78507c and the patches; rerun the identity there, two seeds and all four pairs, any difference
or missing pair stopping the run; then calibrate: C0's long run on roaches and crimtab widened to 3 and on
mite as fitted, at a fixed length (500 sweeps of burn-in, 2,500 kept, 4 chains, three at once), against
the laptop's times for the same runs (95, 92 and 77 seconds, two at once). Their median ratio rescales the
sizing table before the queue starts and converts the box's measured seconds to the table's for the stop.
Then the plan writes look 1's queue, and four workers split it.

**Part B in the harness.** The main run, its scoring and the guard's jobs, written to the plan; the calls
it makes where the plan leaves one open, for the pre-registration:

- *Cells.* 73 cells and 55 folds; the harness writes the cell table (band, rows, S and S_id, identified
  share, mean count, shape, G's columns, R's column and rows, replications) for the pre-registration. BCI
  enters no fitted cell, as re-priced after stage 2. The five shape-below-1 cells take their band's count.
  D4 at S 2 is a small-sample cell, as the design table names it; by the band rule alone its S_id of 1.80
  (over the 87 percent of its rows identified) would put it in the wide band. D6 falls in the middle
  band (1.48). 540 replications at look 1: 485 on truths and 55 folds.
- *Replications.* Narrow 5, middle 5, wide 5 (crimtab's cells 10), small-sample 10, cut from 20 after the
  review: at stage 2's paired spread the band's standard error is about 0.0014 and its noisiest cell's
  0.0033, under the tie margin.
- *Designed cells* not run in stage 2, as their table states: D1's exposure log-uniform over 1 to 1,000 at
  a rate of 0.05; D7's treatment drawn with probability logistic((x1 - x2) / 2), its log rate ratio
  0.3 + 0.7 logistic(1.5 x5), and the rest scaled so that S is 1.
- *G*: per replication, the levels of the first covariate with 2 to 10 distinct values and the quartiles
  of the first numeric covariate with more; each group's log mean count. The columns this picks: roaches,
  treatment and roach1's quartiles; solder, Opening; aids, quarter and delay's quartiles; mite, Substrate
  and SubsDens's quartiles; Arabidopsis, reg; grouseticks, YEAR and HEIGHT's quartiles; epil, trt and
  base's quartiles; Insurance, District; crimtab, finger's quartiles; Traffic, limit and day's quartiles;
  quakes, lat's quartiles; the designed cells, b and x1's quartiles.
- *R*: up to 200 rows a replication, each row's log rate ratio between the two values of a binary column
  and the sample-average rate ratio over them, the two weighted equally; on the designed cells (b) and on
  the six real datasets with a binary covariate: roaches (treatment), epil (trt), Traffic (limit),
  Arabidopsis (amd), solder (Solder) and mite (Topo). Their truths at either level: the additive truth's
  GAM refitted by the export's rule, which reproduces the stored truth exactly, and predicted with the
  column flipped; a widened truth's widening, an affine map of the log mean, applied to that prediction
  too; the BART-fitted truth's fit rerun at its seed with the flipped rows as test rows, which reproduces
  the stored truth exactly since prediction draws nothing. R's cells by band, the shape-below-1 cells
  apart: narrow 4 of 8 (Traffic's two, D5, D8); middle 12 of 14; wide 15 of 26; small-sample 10 of 20
  (roaches' and epil's eight, D4 at both spreads).
- *Readings* are taken from all of a run's draws, as the correction was measured; 500 draws a run are
  stored in single precision for later reweighting (C6 crossed with step 3's choice). A run keeps its
  draws as the chunks it drew, not appended by copy, and is read a block of rows at a time; C0's draws are
  freed before the own runs. Measured on mite with R's rows (2,850 rows) at 5,000 and 10,000 sweeps
  a chain: R's heap peaks at 321 and 413 MB and the process at 850 and 992 MB, growing 1.06 and 1.6
  times as fast as the draws held. Near the cap (100,000 sweeps, 1.8 GB of draws) that is about 2.7 to
  3.7 GB a worker, at most about 15 GB on four workers, against the box's 24 GiB with swap off.
- *Expected coverage*: the long reading against 0.9 N / (N + 1); the default reading against 0.9, since
  what users get includes its Monte Carlo error.
- *The coverage correction, paired.* Steps 2 and 3 read a candidate's gap less C0's, so the correction's
  error that matters is the paired one: a candidate's error less C0's on the same run, rows and share of
  draws, from the correction's runs. At 100 or more effective draws it is +0.0002 (standard error 0.0002),
  every candidate within 0.0004 (chi(4, 1.03) +0.0002, chi(5, 0.9) +0.0004, C4 at 0.6 +0.0001 and at 0.75
  +0.0003, C5 -0.0002). That is the evidence for the own-run threshold of 100 (the second cut of the
  [Amendment after stage 2](#amendment-after-stage-2)); the absolute error's 0.002 holds at the point
  estimates only (above 200 effective draws -0.0015 to -0.0018, about two standard errors). The correction
  was measured on F's training rows; G and R carry it by extension, and P is scored against 0.9.
- *P*: the randomized PIT's uniform for each held-out row is drawn by the cell's, replication's and row's
  own seed, the same for every candidate and both readings.
- *Bounds.* Band-averaged statistics take the paired bootstrap over whole replications within a cell, each
  cell resampled by its own seed, so that every candidate, class and statistic resamples the same
  replications. Cell-level bounds take t on m - 1 degrees of freedom with sd / sqrt(m) over the cell's m
  replications: a cell is worse by the bar when its gap difference less 0.02 exceeds that bound, the
  per-replication difference taken as each reading's miss signed by its mean miss (the gap's linear
  part); the width guard's ratio has the delta method's standard error, sd(IS_c - ratio IS_0) / (sqrt(m)
  mean(IS_0)). The log score per dataset over its 5 folds, C0's less the candidate's, averaged over the 11
  datasets, at its 95 percent bound. P on the simulated cells' held-out rows and on the folds, each
  dataset's folds in its as-fitted band.

**The decision, group-sequential**, at each look on every replication run so far:

- *Boundaries.* Two-sided at 0.05 by Lan-DeMets spending of O'Brien-Fleming type on each side, at the
  realized information fractions, the final look's 1. At equal thirds this gives 3.71, 2.51 and 1.99;
  the looks are never equal here, since the small-sample band and crimtab's cells start at 10
  replications and the rest at 5, so the spending's boundaries replace the classical 3.47, 2.45 and 2.00.
  A look's information fraction is the variance of a band-weighted average over the truth cells at the
  planned maximum (each cell's look-1 replications and two blocks of 5) over its variance at that look,
  each cell's variance a replication estimated at the look, averaged over the classes and fitted
  candidates.
- *Step 2*, per fitted candidate on the long reading: each class's gap less C0's against 0.02, the widest
  cells' gap against 0.02 in each class, and each class's interval score ratio against 1.05 (F, G, R),
  each z = (bar - estimate) / se: passing at the boundary, failing at minus the boundary, pending between.
  No cell worse by the bar, none of the shape-below-1 cells worse, and the log score, as bounded above;
  each narrow cell's width ratio passing at most 1.10 at its upper bound, failing above it at its lower
  bound, pending between. Step 2 fails if any part fails, passes if every part passes, and is pending
  otherwise.
- *Step 3*, per pair of candidates not failing step 2: the lead in rank (the mean of F's and R's gaps less
  C0's) less 0.005 against the boundary: untied at it, tied at minus it, pending between.
- *A look decides* when the recommendation is the same under every way its pending results could resolve
  (a step 2 passing or failing, a pair tied or untied). At the final look pending reads as a fail for
  step 2 and as tied for step 3.
- *Confirmation blocks.* A look that does not decide funds the next look's block: per band, the share of
  each pending statistic's variance a block of 5 replications on that band would remove, summed over the
  pending statistics, per core-hour of the block (the band's measured cost a replication); every band
  within a quarter of the best band's rate, best first, while the stop allows, the last band funded cut to
  its cheapest cells that fit. The folds never grow. A look that can fund no block is final.
- *The stop*, 110 core-hours in the sizing table's units, the box's seconds converted by its calibration.
  At the sizing below, a block on every band costs about 37 and both blocks on every band come to
  about 122: the stop cuts the second wide block, as pre-registered.
- *Step 5*, computed once a look decides: step 3's order (its choice, then its choice among the rest, a
  pending result read as at the final look); the first candidate in it whose default reading's gap less
  C0's is below 0.02 at its 95 percent bound in every class, read once on every replication run; a
  candidate whose default readings are missing waits for them, the decision writing its guard queue. If
  none, the default is kept.
- *Steps 1 and 4.* No candidate passes A with margin and the default fails A (stage 1), so the first
  tie-break and step 4 never apply.

**The widest cells' guard** reads crimtab's cells and the cells widened to S_id 3, whose S runs past A1's
line of 4 on five datasets, all of them at S_id 3 (their spread is over rows whose mean the data barely
identify); the caps hold on all (expected totals at most 19,897, the largest mean at most 9,960):

| cell | S | S_id | cell | S | S_id |
|---|---|---|---|---|---|
| mite widened to 3 | 9.54 | 3 | grouseticks widened to 3 | 3.61 | 3 |
| solder widened to 3 | 5.02 | 3 | crimtab-bart | 3.59 | 2.33 |
| crimtab widened to 3 | 5.01 | 3 | roaches widened to 3 | 3.00 | 3 |
| epil widened to 3 | 4.43 | 3 | Insurance widened to 3 | 3.00 | 3 |
| Arabidopsis widened to 3 | 4.17 | 3 | Traffic widened to 3 | 3.00 | 3 |
| quakes widened to 3 | 3.92 | 3 | crimtab | 2.96 | 2.19 |
| aids widened to 3 | 3.65 | 3 | | | |

**Own runs a replication**, reconciled: stage 2's runs read at C0's stop and the correction's runs, which
ran to it, measure the same rule on different runs. Stage 2's give, a replication, narrow 2.40 (15 runs), middle 0 (5), wide 1.33 (6) and
small-sample 2.05 (20); the correction's 2.25 (4), 0.20 (5), 1.17 (6) and 1.88 (8). The earlier sizing's
wide 1.6 had left mite and solder out of the band, and the amendment's figures are the correction's alone.
Pooled: 2.37, 0.10, 1.25 and 2.00. The sizing uses the pooled rates.

**The smoke**, rerun after these changes on its own seeds: its replications are numbered from 900 and
the main run's go to 30 at most, so no seed is shared; at a tenth of every length (100 sweeps of
burn-in, chunks of 250, a target of 10 effective draws and an own run below 10); scoring off. Look 1 ran
every cell at two replications (the folds at one), each its long and C0's default reading. Each look's
decision ran blind, printing only shapes: every statistic, cell bound, information fraction and boundary,
the enumeration of outcomes, and on a smoke the block allocation whatever the outcome, which funded the
small-sample and narrow bands (28 replications) for looks 2 and 3. Step 5 wrote the guard queue for the
first candidate in its order, waited for the 28 replications each block added, read it once complete,
and moved to the next candidate's queue when that one did not hold. Guard readings of every kind ran:
chi(4, 1.03) on today's sampler on every replication, C4 at 0.6 and 0.75 and C5 on the scratch build
with the move off. 775 jobs on the laptop's two cores, all clean, 1.6 core-hours; 194 own runs (81 chi,
55 floor, 58 half-normal), and no floor run drew a k below its floor; no coverage missing.

**Sizing**, per component, at the amended sizes: stage 2's sweeps a replication by band
(the widened-to-3 wide cells at the mean of the two measured, crimtab's and roaches', 11,250 sweeps a
chain), the move build's seconds a sweep fitted on rows, total count and rows times shape (R squared
0.97), R's test rows at 0.14 of a row each from the smoke, the pooled own-run rates with own runs priced
at C0's length, the readings timed on the smoke, and today's sampler at 1/1.2 of the move build's cost a
sweep for the default readings (stage 2 measured about 1/1.3). Laptop core-hours throughout; the box's
calibration rescales them.

| component | count | core-hours |
|---|---|---|
| C0's long runs | 485 replications and 55 folds | 18.0 |
| own runs | narrow 2.37, middle 0.10, wide 1.25, small-sample 2.00 a replication | 23.3 |
| readings from all draws | | 2.6 |
| C0's default readings | 540 | 2.2 |
| the guard, one candidate's default readings | 540 | 2.2 |
| expected before confirmation | | about 48 |
| a confirmation block of 5 replications | narrow 4.7, middle 3.2, wide 25.2, small-sample 3.8 | 37 a block on every band |
| two blocks on every band | | about 122, past the stop of 110 |
| the smoke, laptop, done | 775 jobs | 1.6 |

By band, before confirmation: narrow 8 cells at 5, 4.7; middle 14 at 5, 2.7; wide 26 at 5 (crimtab's 10),
26.2; small-sample 20 at 10, 7.5; the shape-below-1 cells, 1.3; the folds, 3.6. On the box's four cores the
main run is about 12 hours of wall time at the laptop's speed; its stored outputs about 3 GB. Not
measured: C4 at its floor and C5 on null data have not run at full length, and own runs are priced at
C0's length.

## Amendment before the fresh block

Made 2026-11-02 under dec-B556, after the first look and before any fit on the replications below. Where it
differs from the text above, it governs. The first look's outcome as registered is unchanged.

**The search.** It read only the prior check and stored draws, and no fit was made for it. Screened: df 2 to 4
by 0.25 against scale 1.25 to 3.5 by 0.25, 90 settings, on A1 (the stored forests) and A3 (crimtab's spread,
160,000 forests, standard error about 0.03). Eight pass both. Each was then read from the first look's 500
stored draws a run, by readRun, scoreReading and lookStatistics at the first look's boundary, with no own runs,
so this is a screen and not the rule's reading.

| prior | A1, limit 1 | A3, floor 1 | step 2 | rank | widest width guard | cells worse | narrow F interval score | smallest Kish, largest k-hat |
|---|---|---|---|---|---|---|---|---|
| chi(2.25, 3) | 0.87 | 1.10 | pending | -0.002 | 1.22 Traffic-bart | | 0.982 | 0.01, 0.84 |
| chi(2.5, 2.25) | 0.99 | 1.33 | fail | -0.001 | 1.07 Traffic-bart | Traffic-bart | 0.986 | 0.48, 0.47 |
| chi(2.5, 2.5) | 0.77 | 1.01 | fail | -0.001 | 1.11 Traffic-bart | Traffic-bart | 0.981 | 0.07, 0.64 |
| chi(2.75, 2) | 0.82 | 1.20 | pending | -0.001 | 1.04 Traffic-bart | | 0.994 | 0.52, 0.30 |
| chi(3, 1.75) | 0.75 | 1.13 | pending | -0.001 | 1.02 D5 | | 1.005 | 0.54, 0.31 |
| chi(3.25, 1.5) | 0.79 | 1.19 | fail | -0.002 | 1.07 D5 | epil-w3 | 1.022 | 0.46, 0.31 |
| chi(3.5, 1.25) | 0.98 | 1.66 | fail | -0.000 | 1.14 D5 | epil-n30-3, epil-w3 | 1.049 | 0.22, 0.36 |
| chi(3.75, 1.25) | 0.69 | 1.18 | fail | -0.001 | 1.14 D5 | epil-n30-3, epil-w3 | 1.047 | 0.23, 0.38 |

Rank is the band-weighted gap over F and R less the default's: all eight are tied by step 3's 0.005.

**The pick: chi(2.75, 2)**, by the rule's order. Along the passing settings a larger scale loses on Traffic-bart
and a smaller one on D5's width and epil-w3. Traffic-bart's flag is not noise from the weights: F's gap there rises
with the scale, at -0.001, 0.018, 0.035 and 0.051 from chi(3, 1.75) to chi(2.5, 2.5), and every one of its five
replications has the larger interval score. chi(2.25, 3) escapes the flag only by its noise: its gap is 0.043 and
its weights keep as few as 4 effective draws, so stored draws cannot read it and the rule would ask for own runs
on 101 replications. Its median k is 3.8, twice the default's. That leaves chi(2.75, 2) and chi(3, 1.75), tied, and
the narrow band's interval score takes chi(2.75, 2). No setting between them can be separated from these two on
stored draws, so the search stops. The stored read has overstated epil's losses before, so chi(3.25, 1.5)'s flag may
not stand, but it also loses that tie-break. Not read: any setting off the grid, and own runs for any setting but
these two, as dec-B556's rule has it. chi(2.75, 2) sits close to the cell bar on Traffic-bart's F, at 0.018 against
0.02, which the block retests.

**The block.** Five new replications on each cell that is not a fold, with each cell's replications continuing
from its last, 6 to 10 or 11 to 15, by the harness's seeds; the folds never grow. On each: C0's long run; chi(2.75,
2)'s own run wherever the rule asks for one (Kish below 0.3, k-hat above 0.7, or under 100 effective draws), and on
every narrow replication, since its weight against the default grows as a power of k; C0's and chi(2.75, 2)'s
default-length runs. No other candidate runs: the five stay failed. The harness's code does not change.

**The stop** rises by 27 units to 137 (dec-B556), against 109.9 spent, so 27.1 units, about 19.1 box core-hours.
The block as priced comes to about 29: long and own runs 26.9, default-length runs 2. So mite-bart, 3.1 units, is
cut now, the deviation from "cut to its cheapest cells" being that epil-w3, the costliest cell and the reason the wide
band runs, is kept. Order: default-length runs held in reserve; narrow, middle, small-sample; then the wide band,
epil-w3 first and the rest cheapest first, so mite-w3 comes last. Whatever the stop does not reach stays out and is
named. Runs stop early if a result is settled that no remaining cell could change.

**Two readings, both stated now.** Neither is a final-look reading of the kind the first look's added section used.

1. The rule as registered: the first look's and the block's replications together, the first look's readings of
   chi(2.75, 2) as composed, at the Lan-DeMets boundary for the realized information fraction, about 2.45. This look
   is the last, but its boundary stays the spending boundary and not the 1.96 of a look made final by the stop.
   A statistic still pending reads as not passed.
2. The block alone, free of how the prior was chosen: the same statistics, bars, cell rule and width guard, on the
   new replications, with its z at 1.96. The held-out log score rests on the folds, which do not grow, so it is
   not part of this reading.

Step 5 reads the default-length runs on all replications, and on the block alone. The rule recommends chi(2.75, 2) only
if both readings and step 5 pass. If either fails, the rule keeps the default, and how it failed is reported. Each
reading's numbers are reported either way, and the result goes to the maintainer (dec-B432); nothing is built until
the maintainer rules.

## The defect in R's truth, and its rescoring

Found by both referees after the fresh block (dec-B558). In [`designedBase`](../../benchmarks/R/nbinom-default-study.R)
the designed cells' log mean before scaling was `g(b)`, evaluated inside `with(x, ...)`; `x` has a column `b`, which
masked the argument, so `truthAt(1)` and `truthAt(2)` both returned the truth at each row's observed `b`. The truth
of R's per-row log rate ratio, and of its sample-average rate ratio, was therefore exactly 0 on D1, D2, D3, D4-2,
D4-3, D6 and D8. The generated counts were always right: they are drawn from `truth`, which evaluates `g` at the
observed `b` and so was never affected, and the same holds for F, G and P. Only the stored truth that R is scored
against was wrong. D5, the null cell, has no effect to carry, and D7 builds its truth by its own code; both were right.

The fix keeps `b` outside `with(x, ...)`. A check sourcing the old and the fixed definitions alone, no fit, finds
`truth` and the covariates identical on all nine designed cells; `truthAt(1)` equal to `truthAt(2)` before the fix
on the seven, and different after it, the per-row effect's sd 0.37 (D1), 0.42 (D2), 0.40 (D3), 0.92 (D4-2), 1.39
(D4-3), 0.61 (D6) and 0.12 (D8); the row's observed `b` picking out `truth` again; and D5 and D7 unchanged.

## How the result goes to the maintainer

The finding is written to docs/design/nbinom-default-study.md (Status FINDING): Part A's table and verdicts
under the ruled criterion; Part B's readings by class, band and cell, long and default, with the unsettled
share and the rows below 0.01 apart; the long reading's check; what the rule recommends and why. A last
section, labeled sensitivity, shows what the open call's other criteria would have recommended; it is not an
option offered after the fact. An opus review reads the finding against the pre-registration before it
goes out.

The maintainer is then shown one decision: keep chi(1.5, 2), or change the negative binomial default to the
recommended candidate. Background first, including that the change buys a coherent prior predictive rather
than better intervals; then what each option does to users (every default nbinom fit's draws move; the
help's statement of the default; a family-specific default beside probit's and logistic's, which the
maintainer has said is no cost in itself, dec-B431); the code change and the gates it takes
(posterior-changing: the nbinom exact gates, equivalence scenario and SBC arm, a design note, and for C4 or
C5 an engine slice with its own plan); and the recommendation. A shape prior change (C6) is a separate
decision, shown after. The move's negative binomial half, which waits on this study (collapsed-scale-move's
open call 4), is reported with its measured mixing, for that plan. Nothing is built until the maintainer
rules; the ruling is recorded as a ledger entry quoting the maintainer, and TODO nbinom-default-prior is
updated or closed.

## Open calls for the maintainer

1. The criterion: ruled, dec-B502 ([Amendment after stage 1](#amendment-after-stage-1)).
2. The long reading: recommend long runs on a scratch build of the collapsed move's negative binomial half,
   reweighted across priors (exact; about 65 core-hours against 97 for the grid, an approximation, and 181
   for today's sampler).
3. Ties: recommend a candidate passing A with margin first, then the smaller code change, then narrower
   intervals on narrow data.
4. Weighting: recommend the four spread bands equally (narrow, middle, wide, small-sample wide), cells
   equal within a band.
5. Engine-work candidates (C4, C5) eligible before 1.0-0: recommend yes; about 150 and 300 lines, no ABI
   event.
6. The data-scaled spread (C7): recommend scoring it here.
7. The shape's prior (C6): recommend scoring it here with its weights fixed now, decided separately from the
   spread.
8. Compute: set by the orchestrator's calls ([Amendment after stage 1](#amendment-after-stage-1)), about
   62 core-hours, re-priced from stage 2, stop at about 110; the laptop through stage 2, the main run on
   the box after the k study's first tranche ([Amendment after stage 2](#amendment-after-stage-2)).

## Not verified

- That the scratch move build samples the posterior beyond stage 2's seven cells; its code is reviewed
  again at stage 4.
- Run lengths, weight failures and unsettled shares on small samples, BCI and the widened truths; stage 2
  measures them, and the costs above are projections until then.
- The box's per-sweep times until its calibration rescales the sizing: every cost here is the laptop's.
  Own runs are priced at C0's length; C4 at its floor and C5 on null data have not run at full length.
  R's test rows are priced from quick runs.
- That every widened truth fits under the count rule after re-centering, and that every GAM truth converges
  under the fixed rule (the probes fitted it on all twelve).
- That C6's per-draw weights keep enough effective draws in every cell; a cell where they do not reports C6
  as unscored there.
- That C4 and C5 behave as their reweighted scores say once drawn by an engine step; only the guard's scratch
  build tests that, at the default length.
- That the box is free after the k study's first tranche, about 20 hours from 2026-10-29.
