# nbinom-default-study: is there a better default prior for negative binomial fits?

Status: PLANNED 2026-10-29 under dec-B432. Runs after the merge to main and before 1.0-0. Not run: only the
probes of [Probes](#probes) ran, 2026-10-29, on main 8a78507c. The blind critique has not run, and the
criterion is the maintainer's first open call ([Open calls for the maintainer](#open-calls-for-the-maintainer)).

agent: a blind critique of this plan first (opus); then one opus implementer writes the harness, runs it and
writes the finding; one review of the harness and smoke (sonnet, the gate-only tier of the plan README); one
opus review of the main reading before it goes to the maintainer.
rng: NEUTRAL. No package code changes: a benchmark harness and its outputs. A finalist the package cannot
express today is checked at the default length on a scratch build carrying only its k draw, never landed.
window: after the merge to main, before 1.0-0 (dec-B432). A change of default it supports goes back to the
maintainer; building it is a separate, posterior-changing slice.
budget: harness about 1,300 lines of R (benchmarks/R/nbinom-default-study.R), reusing the first study's
prior predictive draw. Compute about 80 core-hours expected, at most about 110, stop at 125 ([Cost](#cost)).

## Summary for the maintainer

You kept the negative binomial default, k drawn under chi(1.5, 2), judged the study behind it flawed and the
default improvable (dec-B431), and ruled that a redesigned study runs before 1.0-0 with its criterion fixed
first and in absolute terms, count data with wide spreads, chains long enough for k to mix, and candidates
beyond the prior on k (dec-B432). This is that study.

**The criterion**, fixed before any run, has two parts.

- What a dataset drawn from the prior should look like, in absolute terms and with no other package's
  default as the yardstick. In at most 1 prior dataset in 100 should the rows' mean counts span more than a
  factor of about 10^8 (a spread of the log mean across rows, S, above 4); at most 1 in 100 should hold a
  count the package refuses (above 10^6); and no real dataset in the study should be rarer than 1 in 1,000
  under the prior. The widest real data found reach S of 2.9 to 3.8 (a table of criminals' finger length
  by height). The default fails the first two, at about 7 in 100 and 5 to 10 in 100.
- What fits must not lose: interval coverage of what users estimate from count fits (the mean count at a
  row, mean counts by group, rate ratios of a binary covariate, a new count), on truths with known values
  built on real covariates, read as distance from nominal on either side, your rule for calibration, with
  the house's bar of 0.02.

A candidate that passes the first part replaces a default that fails it if its coverage is no worse by the
bar; were the default to pass, a candidate would have to be better by the bar. The numbers in this
criterion are judgments you may want to make yourself; that is the first open call, and the study stores
what it needs to be scored again under another choice without running anything again.

**Compared**: the default; k fixed at 1.5 and at 2; k drawn under chi with 3, 4, 5, 7 or 10 degrees of
freedom at the default's median; the default's form at two-thirds the width; the default with its tail cut
at three places (k kept above a floor, so the forest's spread has a ceiling); a half-normal prior on the
forest's spread at two widths; and a shape prior that puts more weight on small shapes. A different leaf
scale is already among them: the scale and k's prior enter the model only through their ratio. All but the
shape prior are scored from one set of fits. The probes' prior draws show, on six designs, that four pass
the first part (chi with 4 or 5 degrees of freedom, the tail cut at its loosest floor, the narrower
half-normal) and the default does not.

**Chains long enough for k to mix.** Under a drawn prior k mixes slowly: its autocorrelation time was 440
to 8,000 sweeps on seven datasets here (median about 1,550), against tens of sweeps for most of what users
read. With k fixed the trees mostly mix in tens of sweeps. So the study fits each dataset at a grid of fixed
k values and integrates k out afterwards, weighting each grid point by how well it fits (path sampling over
k). That gives the fully mixed posterior under every prior on the spread from one set of fits, paired, and
lets any prior be scored later. In the probes it matched long direct runs wherever those converged; on two
datasets the trees mix slowly even at fixed k, and there neither method settled. So it is checked again
before anything depends on it, and a replication where it does not settle is flagged and read apart. What
users get from a default-length run is read from direct runs, for the default and the chosen candidate.

**Data**: twelve real count datasets with spreads from about 0.1 to about 3.8 (named below), with truths
fitted to each as found and widened to a spread of 3; designed cells for exposure offsets over three orders
of magnitude, mostly zeros, high counts, small samples, no signal and 5,000 rows; and the real outcomes
themselves, scored on held-out counts.

**Stages**: harness and checks; the prior predictive of every candidate (about an hour); the check of the
grid against direct runs; a smoke; the main run on every cell; the default-length check of the chosen
candidate; a confirmation on fresh replications only if the call is close; the finding to you. About 80
core-hours expected, stop at 125. The laptop runs everything before the smoke on two cores while the box
finishes the k study; the box runs the smoke and the main run.

## Goal

After this study each candidate default for `family = "nbinom"` carries a pass or fail on each part of a
criterion fixed beforehand, with the numbers behind it, and either one candidate is recommended to the
maintainer as the new default or the default is recommended kept, with what keeping it costs. The stored
outputs let the maintainer, or a later reader, apply another criterion or another prior on the spread
without new fits.

## What is known

- **The prior as built.** A row's log mean is c + f(x) + o: c is the data's log rate, a transform the fit
  sets (the log of the total count over the total exposure), o the offset, f the forest
  ([nbinom-log-mean.md](nbinom-log-mean.md), [Leaf prior](nbinom-log-mean.md#leaf-prior)). Given k the
  forest's sd at every x is s = A / k, A = 3 ([`defaultLeafScale`](../../R/model.R)), whatever the tree
  count. k is drawn under chi(1.5, 2) ([`resolveLeafHyperprior`](../../R/model.R),
  [`drawsLeafKByDefault`](../../R/spec.R), the conjugate draw [`ChiKHyperprior`](../../src/bartcore/model.hpp)),
  the same prior as sd = invchi(1.5, 1.5): median k 1.91, median s 1.57. The shape r is drawn on the grid
  1 to 50 under gamma(2, 0.1) weights ([`NBShapePrior`](../../src/bartcore/model.hpp)): 15 percent of its
  mass at r of 3 or less, 46 percent at 10 or more.
- **The refused count.** A count above 10^6 is refused at creation and at every response swap
  ([`maximumCount`](../../src/bartcore/model.hpp), dec-B202): the exact Polya-Gamma draw costs one draw per
  trial, about 160 ns, every row and sweep, and the shape's count histogram is sized by the largest count.
- **The first study** ([nbinom-prior-predictive.md](../design/nbinom-prior-predictive.md)). Under the default
  the spread of the log mean across a dataset's rows, S, has median 0.93, 90th percentile 3.1, 99th 15;
  at a mean count of 10, 5 to 10 percent of prior datasets hold a refused count. Its twelve real datasets
  sit at S of 0.12 to 1.4. No candidate fitted better beyond fold noise; k fixed at 2 covered 0.65 against
  0.92 on wide data at n = 30. Its flaws are listed in
  [Limits of this study](../design/nbinom-prior-predictive.md#limits-of-this-study) and each is answered
  here: the criterion is absolute and fixed first ([The criterion](#the-criterion)); the wrongly counted
  generator limit is gone (the cap is the package's, 10^6); the candidates and the rule that picks among
  them are fixed before any fit ([Candidates](#candidates), [The decision rule](#the-decision-rule)); the
  data are wide and include offsets and correlated columns ([Data](#data)); k is integrated out
  ([Chain lengths, and k integrated out](#chain-lengths-and-k-integrated-out)); and the leaf scale and the
  shape's prior are varied.
- **k's mixing under nbinom** ([k-mixing-pg-families.md](../design/k-mixing-pg-families.md)): on datasets
  drawn from the default prior at n = 150, k's autocorrelation time has median 191 sweeps and 99th
  percentile 3,253; its slow cases are mostly-zero data and small true k. On real and wide data it is
  longer ([Probes](#probes)).
- **Cost grows with the counts.** A sweep costs about 0.16 microseconds per unit of y + r summed over rows,
  r the drawn shape, plus about 0.9 microseconds per row, at 75 trees ([Probes](#probes)). Wide spreads make
  large counts: a simulated design with S = 2 and a median count of 5 at n = 1,000 summed to 346,000 (55 ms
  a sweep), and at S = 3 with a right-skewed log mean a count above 10^6 appeared even at a median count of
  0.3. Every wide cell here is sized by its count total.

## The criterion

Fixed before any run. Two parts: A, the prior predictive, decided on prior draws alone in stage 1; B, the
posterior, decided on fits. Both are judgments of the maintainer's; this section states the recommendation,
[Open calls for the maintainer](#open-calls-for-the-maintainer) its alternatives, and
[Stored for re-scoring](#stored-for-re-scoring) what lets the study be scored under any of them later.

### Part A: what a prior dataset should look like

Per prior dataset: S, the sd across rows of c + f(x) (the offset excluded: it is data, not prior); the
largest count; the mean count over the center's exp(c). Cells: each real dataset's own design (its
covariates, n, offset and center), and a simulated grid of n 50, 500, 2,000 and 10,000, centers (mean
counts) 0.1, 1, 10 and 100, and three designs (5 correlated normal columns; 20 uniform columns; 5 normal
and 5 binary columns with an exposure offset spread over three orders of magnitude). A candidate passes A
when, in every cell:

- **A1, spread.** P(S > 4) is at most 0.01. S = 4 means the middle 98 percent of a dataset's rows have
  mean counts spanning a factor of exp(2 x 2.33 x 4), about 10^8. The widest real count data found have S
  of 2.9 to 3.8, by the fit ([Real datasets](#real-datasets)); a prior that goes past 4 in more than 1
  dataset in 100 spends its mass where no data have been seen.
- **A2, fittable.** P(a count above 10^6) is at most 0.01, in cells with a center of 100 or less. A default
  should not generate, in more than 1 dataset in 100, data the package refuses. Cells with a center above
  100 are reported, not read: a mean count of 1,000 at S = 2 reaches the cap without anything being wrong.
- **A3, room.** The S of every as-fitted truth of [Truths on real covariates](#truths-on-real-covariates)
  is at or below the prior's 99.9th percentile of S at its dataset's design: no real dataset is rarer than
  1 in 1,000 under the prior. A1 limits how often the prior proposes spreads no data show; A3 that it does
  not rule out the ones real data do. The two read different tail levels because at the same level, 1 in
  100, they leave no room: the probes found no candidate that passes A1 with a 99th percentile at crimtab's
  design above that table's spread ([Probes](#probes)). crimtab's GAM truth sets the bar at its design
  (3.8); part of its spread lies in corners of the table with no counts, where the data bound the mean only
  from above, so that bar is the high end of what the data require.

Reported beside, not read: the quantiles of S, the share of datasets with more than 90 percent zeros, the
spread of the mean count around exp(c) (the forest's level, which the data-derived center already fixes),
the count total (cost), and each real dataset's percentile.

The default fails A1 (5 to 8 percent of prior datasets above 4 at the probes' six designs) and, by the first
study, A2 (5 to 10 percent); it passes A3. Stage 1 states each candidate's verdict before any fit.

### Part B: what fits must not lose

Estimand classes, what users report from a count regression:

- **F**, the log mean at each row, training and held-out (equivalently the mean count); pointwise.
- **G**, the mean count averaged over a covariate-defined group (a factor's levels, a numeric covariate's
  quartiles); one scalar per group.
- **R**, rate ratios of a binary covariate where the design has one: the per-row log rate ratio
  f(x, 1) - f(x, 0) and the sample-average rate ratio, the ratio of average mean counts at 1 and at 0.
- **P**, a new count at a held-out row: 90 percent predictive coverage by randomized PIT,
  u = F(y - 1) + V (F(y) - F(y - 1)).
- The shape r is reported, not read.

Measures, per cell, class and candidate: coverage of 90 percent intervals (50, 80 and 95 read beside) of the
truth; the cell's gap, its distance from nominal on either side corrected for noise as
[Decision rule](binary-interval-study.md#decision-rule) defines it (dec-B455); the interval score; RMSE of
the log mean; and on real outcomes the held-out log predictive density per row.

**The long reading** is each candidate's posterior with k integrated out
([Chain lengths, and k integrated out](#chain-lengths-and-k-integrated-out)): the prior's own effect, with
mixing taken out, as dec-B432 asks. **The default reading** is a direct run at bart()'s defaults (4 chains
of 500 burn-in and 500 kept, 75 trees), what users get; every candidate's per-sweep time is the same to
within a few percent (the k draw is one gamma variate a sweep), so equal sweeps are equal time.

Bars, the house's (September's, reused by the k study,
[The bars](../design/binary-k-fixed-or-drawn.md#the-bars)): 0.02 of coverage gap; 0.01 nats per row of
held-out log score. Bounds: a class average is "worse by the bar" when the upper end of its 95 percent
interval (bootstrap over whole replications, paired) reaches the bar; a single cell when its point
estimate exceeds the bar by 1.96 of its standard errors.

### The decision rule

1. **A**, from stage 1, before any fit: each candidate passes or fails.
2. **Non-inferiority, long reading**: in each of F, G, R and P, the candidate's average gap less the
   default's is below 0.02 at its upper bound; no cell is worse by the bar; the held-out log score on real
   outcomes is not worse by 0.01 nats per row at its upper bound; the interval score is not worse by more
   than 5 percent in any class.
3. **If the default fails A** (expected): among candidates that pass A and step 2, the one with the
   smallest average gap over F and R is recommended; candidates within 0.005 of it count as tied, and a tie
   goes to the smaller code change ([Candidates](#candidates)), then to the heavier tail (the larger 99.9th
   percentile of S pooled over A's cells, the draws nearest the default's). If none passes, the default is
   recommended kept, with what passing A costs in coverage.
4. **If the default passes A**: a candidate is recommended only if it also improves F or R by at least
   0.02 at its lower bound; otherwise the default is kept.
5. **The default reading guard**: the recommended candidate's default reading is not worse than the
   default's by the bar in any class (step 2's bounds). If it is, the next candidate by step 3's order is
   checked. A candidate better than the default only at the default reading is reported as a mixing
   finding for TODO k-mixing-pg-families, not as a better prior.
6. **Close calls are confirmed**: if the recommended candidate's step 2 margin in any class, or its lead
   in step 3 or 4, is under half the bar, a confirmation on fresh replications (a third of the main run's
   count, every cell) repeats steps 2 to 5 on the long and default readings, and decides.

### Stored for re-scoring

Everything the rule reads is recomputable from stored outputs, under any of the alternatives of the first
open call and under any prior on the forest's spread:

- Prior predictive, per cell and draw, for a unit-spread forest: S at unit spread, the mean, the row
  quantiles at 101 levels and the 20 largest values of f. S and the largest count under any prior on the
  spread follow from these by rescaling (S is proportional to the spread), with counts drawn at the top
  rows (checked against full draws on a subset in stage 1).
- Fits, per replication and grid point: per draw the leaf count, leaf sum of squares and shape; per row
  the posterior mean, sd and 19 quantiles (2.5 to 97.5 percent) of F; G and R draws thinned to 200; P's
  cdf at y and y - 1 averaged over draws; the truth; the grid's acceptance. Default readings store the same
  per run, plus k and shape draws, split-Rhat and effective draws.
- About 3 GB in all. A prior on the spread is scored by recomputing grid weights; a prior on the shape
  among the pre-registered ones by weights stored per draw.

## Data

### Real datasets

All available offline: base R's datasets and recommended packages, and CRAN packages installed on the
laptop, exported once to .rds files for the box. S is the sd across rows of the fitted log mean, offset
excluded, from a main-effects NB GLM (MASS's glm.nb), an NB GAM (mgcv's nb() family) and, where measured,
a dbarts fit at k fixed at 1; theta is the GAM's fitted shape ([Probes](#probes)).

| dataset (package) | n | total count | zeros | S: GLM, GAM, BART | theta | note |
|---|---|---|---|---|---|---|
| roaches (rstanarm) | 262 | 6,720 | 0.36 | 1.32, 1.37, - | 0.33 | exposure offset |
| solder (survival) | 900 | 4,977 | 0.32 | 1.36, 1.47, 1.66 | 5.2 | |
| aids (boot, dud == 0) | 465 | 6,215 | 0.21 | 1.57, 1.71, 1.72 | 4.5 | |
| mite, long form (vegan) | 2,450 | 9,800 | 0.57 | 1.55, 3.31, 1.97 | 0.73 | 35 species by 70 cores |
| Arabidopsis (lme4) | 625 | 18,727 | 0.20 | 1.07, 1.08, - | 0.60 | |
| grouseticks (lme4) | 403 | 2,567 | 0.31 | 1.36, 1.37, - | 0.77 | |
| epil (MASS) | 236 | 1,950 | 0.10 | 0.71, 0.81, - | 2.8 | |
| Insurance (MASS) | 64 | 3,151 | 0.02 | 0.31, 0.31, - | above 10^5 | holders offset; Poisson |
| crimtab, long form (datasets) | 924 | 3,000 | 0.67 | 2.92, 3.80, 3.42 | 72 (GLM 0.78) | finger length by height |
| BCI, long form (vegan) | 11,250 | 21,457 | 0.60 | 1.79, -, - | 1.6 (GLM) | 225 species by 50 plots |
| Traffic (MASS) | 184 | 3,965 | 0.00 | 0.12, 0.15, - | 11 | |
| quakes (datasets) | 1,000 | 33,418 | 0.00 | 0.49, 0.49, - | 20 | stations reporting |

Stage 0 fits the missing BART spreads and BCI's GAM (the GLM's formula with a smooth of location). Four
have a shape below 1, the shape grid's floor (TODO approximate-distributions); their cells are also read as
a split, since there the model cannot reach the truth's shape.

Excluded from fits, and reported against each candidate's percentiles: Migration (carData; S 1.55 to 1.57;
its counts total 830,000, about 1.5 hours a replication at the measured cost).

### Truths on real covariates

For each real dataset of the table, two truths, fixed once from the real outcome and stored:

- **Additive**: the log mean of an NB GAM (mgcv's nb() family; smooths of numeric covariates, factors as
  given, the formula fixed per dataset in the harness), with its fitted theta as the true shape.
- **BART-fitted**: the posterior mean log mean of a long dbarts fit at k fixed at 1 (s = 3, wider than every
  candidate's center, so the truth is not one candidate's own shape), with the GAM's theta.

Each at two spreads: **as fitted**, and **widened** to S = 3, about the widest real data's, by scaling the
truth's departure from its mean; a truth already at 3 or more (crimtab's two, mite's GAM) has no widened
cell. Widened truths whose expected count total exceeds 50,000 or whose expected largest count exceeds 10^5
have their center lowered until they do not, and carry the resulting mean count as a label; one
beyond-edge truth per dataset at S = 4.5, past A1's 4, is run on two replications and reported, not read.
Counts are drawn as rnbinom(size = theta, mu = exp(truth + offset)), non-integer theta included: users'
shapes are not on the grid. Held-out rows are a fixed fifth of each dataset's rows.

### Designed cells

Simulated covariates, each isolating one feature; truths are smooth functions with interactions, scaled to
the stated S:

| cell | n | feature | S | center | shape |
|---|---|---|---|---|---|
| D1 offset | 500 | exposure log-uniform over 1 to 1,000 | 1 | rate 0.05 | 3 |
| D2 zeros | 1,000 | mostly zeros | 1 | 0.2 | 2 |
| D3 high | 300 | high counts | 1 | 100 | 10 |
| D4 small | 60 | small sample | 2 | 5 | 3 |
| D5 null | 500 | no signal | 0 | 5 | 5 |
| D6 large | 5,000 | many rows | 1.5 | 2 | 3 |
| D7 effect | 1,000 | binary treatment, heterogeneous log rate ratio 0.3 to 1 | 1 | 5 | 3 |

Every designed cell has 5 correlated normal covariates (correlation 0.5 between neighbors) and one binary
column; D7's treatment depends on them, as an observational design would.

### Real outcomes

Each real dataset's own outcome, scored by 5-fold held-out log predictive density per row and P's coverage
of the held-out counts, two repetitions. No truth is involved, so it enters only P and the log-score guard.

## Candidates

Every candidate below except the shape prior is a prior on the forest's spread s = 3 / k, so all are scored
from the same fixed-k grid ([Chain lengths, and k integrated out](#chain-lengths-and-k-integrated-out)), and
a further spread prior costs nothing to score. So each form is sampled at several widths, placed from the
probes' prior draws to straddle A's bounds. A different leaf scale A is not a separate candidate: s = A / k,
so A and chi's scale enter only through their ratio, and a rescaled default is chi(1.5, c), the family C3
samples. A scale read from the data, as gaussian's response range is, was considered and is not a
candidate: a count response does not separate the spread of the log mean from the overdispersion without
fitting a model, so no statistic of y alone measures what gaussian's range does.

| | candidate | s's prior | why | code change if chosen |
|---|---|---|---|---|
| C0 | default | k ~ chi(1.5, 2), s = 1.5 / chi(1.5), median 1.57 | | none |
| C1 | k = 1.5, k = 2 | s fixed at 2 or 1.5 | no tail at all, no k to mix | a family default for k in place of [`drawsLeafKByDefault`](../../R/spec.R)'s flag; R about 30 lines, help, prior-defaults.md, tinytest; fits lose their k channel |
| C2 | chi(3, 1.25), chi(4, 1.03), chi(5, 0.9), chi(7, 0.75), chi(10, 0.62) | the default's median k, 1.9, thinner tails | the first study's best-looking direction, now chosen before fits | as C1 |
| C3 | chi(1.5, 3) | the default's form at two-thirds the width | narrows bulk and tail together | as C1 |
| C4 | chi(1.5, 2) with k at least 0.6, 0.75 or 1 | the default, s capped at 5, 4 or 3 | cuts only the tail A1 objects to; leaves the bulk exactly | [`ChiKHyperprior`](../../src/bartcore/model.hpp) gains a floor (exact by rejection); [`chi`](../../R/model.R) and invchi gain a bound; the header's hyperprior struct a field (an ABI event); about 200 lines and an exact gate |
| C5 | half-normal on s, median 1.3 or 1.57 | thin tail without a hard cap; 1.57 is the default's median | the standard weakly informative form | a non-conjugate k step (slice on log k), a new constructor, the header, docs that say invchi is the only form; about 400 lines and an exact gate |
| C6 | the shape on the grid with log-uniform weights | 30 percent of mass at r of 3 or less (15 under gamma(2, 0.1)) | real counts mostly have small shapes ([Real datasets](#real-datasets)); crossed with the default and with step 3's choice | the weights in [`NBShapePrior`](../../src/bartcore/model.hpp) and its mirrors in benchmarks/R; help; about 40 lines |

Every change of default is posterior-changing for nbinom fits: the exact gates' nbinom entries, the nbinom
equivalence scenario and SBC arm re-derived, and a design note, per the
[RNG classes and their gates](README.md#rng-classes-and-their-gates).

Not candidates, with the reason: the tree prior and tree count (the tail is in s; S's dependence on the
tree prior is reported from stage 1); a real or sub-1 shape (not built, TODO approximate-distributions);
the grid's cap at 50 (Poisson-like data at large means; a separate question); the center c (data-derived;
moving it is a different model, not a prior).

## Chain lengths, and k integrated out

**What direct runs would take.** k's autocorrelation time under the default, measured here
([Probes](#probes)), is 440 to 3,700 sweeps on six datasets of 262 to 1,000 rows, median about 1,550 with
mite's, and about 8,000 on mite (2,450 rows; 7 effective draws, so unreliable). chi(5, 0.9) does not
reliably shorten it (2,450 on roaches against 1,200; 4,600 on mite against 8,000). On most datasets the
estimands mix in tens of sweeps whatever k does, but not on all: on solder their autocorrelation times are
1,100 to 1,900 sweeps and on mite up to 6,100, and a default-length run's 90 percent intervals hold 0.73 to
0.79 and 0.54 to 0.92 of the long run's draws. About 100 effective draws of k take 155,000 sweeps at the
median, 370,000 at solder's and 800,000 at mite's; for the default and one other candidate over the main
run's replications that is about 350 core-hours. Direct long runs are therefore used only to check the grid
(stage 2).

**k integrated out on a grid.** For each replication:

- Fixed-k fits on a grid in log k, as two path chains, one walking the grid upward and one downward, each
  moving k with `setLeafPrior` and burning 100 sweeps after each change before keeping 1,000. At fixed k
  only the trees and the shape move. Per kept draw: the leaf count L and the sum of squared leaves SS, read
  from the saved trees, the shape, and the estimands.
- Path sampling (Gelman and Meng 1998): since only the leaves' prior depends on k,
  d log p(y | k) / d log k is the posterior mean at that k of L - k^2 m SS / A^2, m the tree count and A the
  anchor, whose prior mean is 0. Its grid means, integrated in log k (cubic Hermite, the slopes known at
  every point), give log p(y | k) up to a constant, with a Monte Carlo standard error from the draws'
  effective count.
- Any prior on the spread, pi(k), then has posterior weights proportional to pi(k) p(y | k) on a fine
  sub-grid; each estimand's posterior is the grid's draws mixed by those weights. A fixed-k candidate is its
  grid point. A prior on the shape is applied by per-draw importance weights, its ratio to gamma(2, 0.1),
  which also corrects p(y | k).
- The grid: step 0.1 in log k (the posterior sd of log k measured here is 0.10 to 0.20, so the step is at
  most one sd); placed adaptively, the downward chain starting at k = 4 and walking down until log
  p(y | k) has fallen 8 below its running maximum and every candidate's weight on the last point is under
  0.5 percent, the upward chain retracing it; k = 1.5 and 2 are always on it. In the probes this was 8 to
  16 points; with every candidate's range covered, about 14 on average.
- **Acceptance, per replication.** At every grid point holding at least 1 percent of some candidate's
  posterior weight, the two chains' integrand means agree within three combined standard errors and the
  integrand has at least 50 effective draws. A replication that fails is run once more at those points
  with 5,000 kept sweeps a chain; one still failing is kept, flagged, and read only in a sensitivity
  reading beside the main one. A cell enters the long reading when at least 80 percent of its replications
  are accepted; otherwise its long reading is reported as unsettled and the cell enters the rule through its
  default reading only. The probes predict failures where the trees mix slowly even at fixed k (solder,
  mite, single rows of the wide design: 3 to 50 effective draws of the integrand from 2,000 kept), so the
  share of cells left unsettled is reported as a finding of its own.

The probes' check of this on five datasets is in [Probes](#probes); stage 2 checks it on 16 replications
before anything depends on it.

**The default reading** is a direct run at bart()'s defaults and needs no mixing of k: it is what users
get. **The grid check's direct runs** are 4 chains of 2,000 burn-in and 40,000 kept sweeps, about 100
effective draws of k at the median autocorrelation time, 45 at solder's and 20 at mite's.

## Probes

Run 2026-10-29 on main 8a78507c, on the laptop (arm64 macOS), two cores, about 1.6 core-hours in all; the
scripts and outputs are not checked in.

**Spreads of real data.** S from the three fits of [Real datasets](#real-datasets) on thirteen datasets
(the twelve and Migration), the dbarts fit (2 chains of 1,000 burn-in and 1,000 kept at k = 1) on four.
Spreads reach 2.9 to 3.8 on crimtab and 1.6 to 3.3 on mite; the GAM's are the widest where most counts are
zero, as the fit's mean over an all-zero region is bounded by the data only from above. Four of the
thirteen have a GAM shape below 1; on roaches and mite every kept draw of the shape sat at the grid's
floor of 1.

**The prior's spread by candidate.** 2,000 forests drawn from the prior at unit spread (k = 3) on six
designs (crimtab, roaches, solder, aids, mite, and 5 correlated normals at n = 1,000) give S at unit
spread, median 0.48 to 0.67; S under a prior on the spread is that times the spread, combined here with
400,000 draws of the spread. They reproduce the first study's default (median 0.91, 90th percentile 3.1,
99th 14 at crimtab's design). A1 and A3 at the recommended levels, on these six designs only (A3 against
crimtab's GAM truth, 3.80, and mite's, 3.31):

| candidate | P(S > 4) over the six designs | 99.9th percentile of S, crimtab's design | A1 | A3 |
|---|---|---|---|---|
| C0, chi(1.5, 2) | 0.054 to 0.084 | 67 | fails | passes |
| chi(3, 1.25) | 0.009 to 0.020 | 9.8 | fails | passes |
| chi(4, 1.03) | 0.003 to 0.0085 | 6.4 | passes | passes |
| chi(5, 0.9) | 0.001 to 0.004 | 4.8 | passes | passes |
| chi(7, 0.75) | 0.0002 to 0.0014 | 3.7 | passes | fails |
| chi(10, 0.62) | under 0.0004 | 3.0 | passes | fails |
| chi(1.5, 3) | 0.029 to 0.046 | 45 | fails | passes |
| C0 with k at least 0.6 | 0.0006 to 0.0067 | 4.5 | passes | passes |
| C0 with k at least 0.75 | under 0.0024 | 3.7 | passes | fails |
| C0 with k at least 1 | under 0.0004 | 2.9 | passes | fails |
| half-normal, median 1.57 | 0.004 to 0.019 | 5.8 | fails | passes |
| half-normal, median 1.3 | 0.001 to 0.0074 | 4.9 | passes | passes |
| k fixed at 1.5, at 2 | 0 | 2.5, 1.9 | passes | fails |

A1's threshold of 4 was set from the real spreads before this table; A3's level of 1 in 1,000 after it,
since at 1 in 100 no candidate passes both. Part A is a property of the prior alone, so the maintainer
fixing it sees this table; what must be blind is B, and no fit's coverage has been read.

**Mixing, direct runs.** Seven datasets: roaches, solder, aids, mite, and three simulated with the 5
correlated normals and a binary column at n = 1,000 ("wide": S = 2, median count 0.3; "narrow": S = 0.3,
5; "zeros": S = 1, 0.3). Each 4 chains of 16,000 sweeps (1,000 burn-in) under the default, three also
under chi(5, 0.9) and at k fixed at 2. Estimands: the mean log mean over rows, the rate ratio where there
is a binary column, and the log mean at the first row and at the row with the largest count. "Default
length" is the share of the long run's draws inside the 90 percent interval from sweeps 501 to 1,000 of
each chain, as bart() would report it.

| dataset | n | total count | ms a sweep | k's autocorrelation time | slowest estimand's | default length |
|---|---|---|---|---|---|---|
| roaches | 262 | 6,720 | 1.5 | 1,200 | 39 | 0.85 to 0.89 |
| aids | 465 | 6,215 | 3.0 | 1,140 | 38 | 0.88 to 0.91 |
| narrow | 1,000 | 5,297 | 2.8 | 440 | 36 | 0.87 to 0.90 |
| zeros | 1,000 | 692 | 2.8 | 1,630 | 95 | 0.88 to 0.92 |
| wide | 1,000 | 18,792 | 5.2 | 1,550 | 590 | 0.69 to 0.92 |
| solder | 900 | 4,977 | 8.1 | 3,720 | 1,910 | 0.73 to 0.79 |
| mite | 2,450 | 9,800 | 3.6 | 8,000 (unreliable) | 6,110 | 0.54 to 0.92 |

Under chi(5, 0.9), k's autocorrelation time was 2,450 on roaches, 1,330 on wide and 4,600 on mite. At k
fixed at 2 the estimands' fell to 5 to 30 sweeps on roaches and zeros, and to 220 at wide's first row (590
drawn).

**The grid.** On five datasets: two path chains, step 0.1 in log k, 100 sweeps of burn-in after each change
and 1,000 kept a chain per point, 8 to 16 points, 55 to 210 seconds each. The units check passed: under the
prior at fixed k the integrand's mean was -0.50, -0.39 and -0.88, standard errors about 0.5, at k 0.8, 2
and 4 (1,500 draws each).

| dataset | k's 5th, 50th, 95th percentiles: grid | direct (k's effective draws) | direct draws inside the grid's 90 percent intervals | integrand effective draws a point |
|---|---|---|---|---|
| roaches | 0.95, 1.26, 1.83 | 0.92, 1.31, 1.82 (48) | 0.895 to 0.907; 0.890 to 0.905 under chi(5, 0.9) | 5 to 220 |
| aids | 1.43, 1.75, 2.09 | 1.37, 1.73, 2.12 (55) | 0.890 to 0.898 | 4 to 98 |
| wide | 0.60, 0.71, 0.84 | 0.61, 0.72, 0.86 (39) | 0.898 to 0.908 for the mean and rate ratio; 0.825 to 0.854 for single rows | 3 to 50 |
| solder | 0.96, 1.13, 1.38 | 1.05, 1.28, 1.56 (17) | 0.888 for the mean; 0.844 to 0.893 for single rows | 3 to 44 |
| mite | 0.57, 0.68, 0.78 | 0.52, 0.64, 0.79 (7) | 0.567 to 0.901 | 3 to 8 |

Where the direct run converged (roaches, aids, wide's mean and rate ratio) the grid matches it within its
noise, under both priors. Where the trees mix slowly even at fixed k (solder, mite, wide's single rows),
neither the grid at this length nor 64,000-sweep direct runs settle, and they disagree (solder's median k
1.13 against 1.28). Hence the per-replication acceptance and stage 2. The slow trees are a finding for TODO
tree-mixing-proposals, not for this study's priors.

**Cost.** Milliseconds a sweep, fitted to the seven direct runs, are about 0.16 microseconds times the sum
over rows of y + r, plus 0.9 microseconds a row, at 75 trees; r is the drawn shape, so Poisson-like data cost
through n r: solder's 8 ms is 900 rows at a shape near 46.

**A package quirk met on the way**: `dbarts()` given `offset.test = v` with `v` an existing variable whose
value is NULL errors "cannot find test offset 'v'"; the probes omitted the argument instead. Reported, not
part of this study.

## Cost

Per-sweep times from the probes' cost model; the box's are assumed equal (Not verified). A replication: the
grid, about 14 points x 2 chains x 1,100 sweeps, 30,800; the acceptance's reruns, assumed for a third of
replications at six points, about 16,000 on average; C0's default reading, 4,000. About 51,000 sweeps.

| stage | where | core-hours |
|---|---|---|
| 0 harness and checks | laptop | about 2 |
| 1 prior predictive, 20,000 prior datasets a cell | laptop | about 1 |
| 2 grid check: 16 replications, the grid and direct runs of C0 and chi(5, 0.9) | laptop | about 6 |
| 3 smoke | box | about 1 |
| 6 main run | box about 60; laptop about 5 | about 65 |
| 7 guard, one candidate (the second only if the first fails) | box | about 5 (10) |
| 8 confirmation, only if called, with its guard | box | about 23 |
| expected: no confirmation, one guard | | about 80 |
| at most | | about 110 |

The main run: each real dataset's truths, as fitted and widened, at 15 replications (BCI 8), and two
beyond-edge replications a truth, about 58; the widened cells dominate, at 50,000 counts about 7 ms a
sweep, and BCI's 11,250 rows cost about 11. The designed cells at 20 replications, about 7. The laptop
takes the cells under 2 ms a sweep (the as-fitted truths of the small datasets, D1, D2, D4, D5) and those
datasets' real-outcome folds.

**Stop**: 125 core-hours over all stages, counted from the harness's own timing logs; at the stop the
reading is taken on what has run and the finding says so. If the smoke's projection exceeds 110, cuts are
made in this order before the pre-registration commit: the beyond-edge truths; BCI to 4 replications; the
widened truths' count cap from 50,000 to 25,000 (lower centers); real-cell replications from 15 to 10.

**Where and when.** Stages 0 to 2, about 9 core-hours, run on the laptop's two cores, about 5 hours of wall
time, inside the roughly 20 hours the box needs for the k study. The smoke and the main run go to the box's
four cores after that, about 16 hours; a confirmation adds about 6.

## Stages and gates

Each gate is a stop: a failure is a harness defect to fix, not a finding, unless it says otherwise.

0. **Harness and checks** (laptop). Data export, truths, the harness, and these checks:
   - **Units**: under the prior at fixed k, the mean of L - k^2 m SS / A^2 (the path integrand, m trees, A
     the anchor) is 0 within three standard errors at three k values; and in a drawn-k run, the posterior
     mean of k^2 equals the mean over draws of its conjugate conditional mean computed from the stored L and
     SS, so the leaves read from saved trees are the ones the engine conditions on.
   - **Path chains agree**: on every smoke replication the upward and downward chains' integrand means
     agree within three combined standard errors at every grid point; otherwise the burn-in after each k
     change is doubled.
   - **A grid point is a fixed-k fit**: the mixture for C1 equals the grid point's own draws.
   - **The prior draw is the fitted model's prior**: cut points, center c and s = A / k as bart() sets them
     (the first study's checks, reused); the top-rows shortcut for the largest count agrees with full draws
     on 1,000 datasets a cell.
   - **Truth-blind**: no grid placement, interval, group definition or acceptance reads the truth.
1. **The prior predictive** (laptop): every candidate on every A cell and real design, 20,000 prior
   datasets a cell; A's verdicts and the real datasets' percentiles written to the harness output and quoted
   in the pre-registration. No fit.
2. **The grid check** (laptop): 16 replications (an as-fitted additive truth on each real dataset but BCI,
   Insurance and Traffic; a widened one on roaches, aids and solder; D2, D4, D5 and D6) also run directly
   under C0 and chi(5, 0.9), 4 chains of 2,000 burn-in and 40,000 kept sweeps. Read on the replications
   whose direct runs converged (split-Rhat of k under 1.05 and at least 50 effective draws of k), at least
   10 of the 16. Passes when, on each of those, k's 5th, 50th and 95th percentiles from the grid are within
   two Monte Carlo standard errors of the direct run's, and, pooled, the share of the direct run's draws
   inside the grid's 90 percent intervals of F, G and R is within 0.90 +- 0.015 with no replication below
   0.85. If it fails, or fewer than 10 converge, the long reading is dropped: every candidate the package
   expresses today (C0 to C3) is read directly at the default length, at about the grid's cost; C4 and C5
   are unscored; the finding says the long reading was not obtained.
3. **Smoke** (box): every cell, one replication, a tenth of every length, every reading and the
   real-outcome folds; the per-cell timing and acceptance rate that the pre-registration's counts and costs
   use.
4. **Review** of the harness and smoke; fixes.
5. **Pre-registration commit** (below), before any study fit.
6. **Main run** (box; the laptop takes cells under 2 ms a sweep on its two cores): every cell's
   replications, grid and C0's default reading; the real-outcome folds. At 10 replications a cell the
   reading is taken once: if every candidate's step 2 verdict clears its bound by at least the bar's half
   and step 3's leader leads the next untied candidate by at least 0.01, the remaining replications are
   skipped (the house's rule to stop once a result is settled).
7. **Reading and guard**: steps 1 to 4 of [The decision rule](#the-decision-rule); then the recommended
   candidate's default reading on the main run's replications (C4 or C5 on a scratch build carrying only its
   k draw), and the next one's only if it fails; step 5.
8. **Confirmation**, only if step 6 of the rule calls it.
9. **Finding**, its opus review, a ledger entry for the calls made running it, and the result to the
   maintainer ([How the result goes to the maintainer](#how-the-result-goes-to-the-maintainer)).

## Pre-registration

Committed to this plan, with the harness at that commit, after stages 0 to 4 and before any stage 6 fit:
the criterion and its numbers; the decision rule; the candidates and their exact priors; the grid (range,
step, extension rule, burn-in and draws per point, acceptance); the replication counts and the truths'
formulas and seeds; stage 1's A verdicts; stage 2's grid-check result; the per-cell costs and the projected
total. Already fixed by this plan and not changed by stages 0 to 4: the criterion, the bars and bounds, the
rule, the candidate list, the cells, the stop and the order of cuts. What stages 0 to 4 may change is
listed, and only that: the grid's step and burn-in (to pass stage 0's agreement check), replication counts
(down, to fit the stop, by the order of cuts), and a widened truth's center (by the stated rule). Any other
change after a fit is a deviation, disclosed in the finding.

## How the result goes to the maintainer

The finding is written to docs/design/nbinom-default-study.md (Status FINDING): A's verdicts with each
candidate's quantiles of S and refused share; B's readings by class, cell and spread, long and default;
the grid check and the share of cells left unsettled; what the rule recommends and why; and what the
alternatives of the first open call would have recommended, computed from the stored outputs. An opus
review reads it against the pre-registration before it goes out.

The maintainer is then shown one decision: keep chi(1.5, 2), or change the negative binomial default to the
recommended candidate. Background first; then what each option does to users (every default nbinom fit's
draws move; the help's statement of the default; a family-specific default beside probit's and logistic's,
which the maintainer has said is no cost in itself, dec-B431); the code change and the gates it takes
(posterior-changing: the nbinom exact gates, equivalence scenario and SBC arm, a design note, and for C4 or
C5 an engine slice with its own plan); and the recommendation. A shape prior change (C6) is a separate
decision, shown after. Nothing is built until the maintainer rules; the ruling is recorded as a ledger
entry quoting the maintainer, and TODO nbinom-default-prior is updated or closed.

## Open calls for the maintainer

1. The criterion ([The criterion](#the-criterion)): recommend A1 to A3 as written (S above 4 at most 1 in
   100; a refused count at most 1 in 100 at centers up to 100; every as-fitted real truth within the
   prior's 99.9th percentile), B on F, G, R and P at the 0.02 bar, read long with the default-length guard,
   and the asymmetric rule (non-inferiority suffices against a default that fails A). Alternatives: A
   reported only, with a candidate needing to be better on B; A3 at 1 in 100, or A1 at S above 3 (by the
   probes, no candidate then passes A); A3 on the GLM's spreads, the narrowest fits (more candidates pass);
   B read on the default length only; 95 percent intervals.
2. Choosing among candidates that pass: recommend the best coverage over F and R, ties within 0.005 to the
   smaller code change and then the heavier tail. Alternatives: the smallest change to the default's draws
   first (C4); the smallest code change first (C1 to C3).
3. Weighting truths: recommend every cell equal within a class, as-fitted and widened truths equally.
   Alternative: weight by how common each kind of count data is in use, which the maintainer would set.
4. The long reading by k integrated out on a grid, checked against direct long runs, with the default-length
   fallback if the check fails: recommend it. Alternative: direct long runs for C0 and chi(5, 0.9), and C1
   at its fixed k, about 350 core-hours more, the other candidates unscored.
5. Candidates needing engine work (C4, C5) eligible before 1.0-0: recommend yes, nothing being frozen
   before release, at their stated cost. Alternative: score them, but recommend only from C1 to C3.
6. The shape's prior (C6): recommend scoring it here and deciding it separately from the spread.
   Alternative: drop it, leaving gamma(2, 0.1) unexamined.
7. Compute: recommend about 80 core-hours, stop at 125, the main run on the box after the k study.

## Not verified

- That path sampling over k matches long direct runs beyond the probes' five datasets, at the large-n and
  mostly-zero cells, and where the trees mix slowly, where in the probes neither method settled; stage 2
  decides.
- That A's verdicts hold beyond the probes' six designs; the probes' 99.9th percentiles rest on 2,000
  unit-spread draws a design. Stage 1 decides.
- That the acceptance's reruns are needed in about a third of replications, as costed; the smoke measures
  it.
- The box's per-sweep times: every cost here is the laptop's, scaled by the probes' per-trial and per-row
  rates.
- That reading leaves from saved trees every 250 draws stays a small share of a sweep's cost at 11,250 rows
  (BCI).
- That every widened truth fits under the count rule after re-centering, and that every GAM truth converges
  with the formula the harness fixes.
- That the shape prior's importance weights (C6) keep enough effective draws in every cell; a cell where
  they do not reports C6 as unscored there.
- That C4 and C5 behave as their reweighted scores say once drawn by an engine step; only the guard's scratch
  build tests that, at the default length.
- The per-replication spreads behind the replication counts.
- That the box is free after the k study, about 20 hours from 2026-10-29.
