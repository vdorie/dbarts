# Samplers for the slow tail of k and the fitted values under probit

Status: SURVEY, 2026-10-09 (dec-B415). Measured on prototypes kept
outside the tree; nothing built. Probit throughout, logistic where a candidate
applies; tree-structure stickiness is out of scope (TODO
tree-mixing-proposals).

## Summary

The probit rescaling step left a tail: on the calibration census k's
autocorrelation time has a 99th percentile near 900 sweeps and the
fitted values' near 1000. The tail has one cause. The step moves k,
the leaves and the latent variables together, and because the latents
come along, how far it can move in one sweep is set by how tightly the
latents pin the fit's scale: about 1/sqrt(2n) in log k, 0.058 at n =
150, whatever the data. On the slow datasets the posterior of log k is
0.7 to 0.9 wide, so the chain needs hundreds of sweeps to cross it.

Only moves that integrate the latents out of the scale update escape
that limit:

- **The collapsed rescaling move** (multiply every leaf by a factor and
  divide k by it, the factor drawn against the probit likelihood with
  the latents integrated out, then redraw the latents). Taken every
  sweep: k's 99th percentile 879 to 124 sweeps, worst 959 to 129; the
  fitted values' 99th 1003 to 64. It costs 70 percent of a sweep at the
  census size and 50 percent at 5000 rows and 200 trees. **Taken every
  fourth sweep it keeps almost all of that** (on the ten slowest
  datasets k's worst 160, the fitted values' worst 45) for 12 to 15
  percent of a sweep. Exact: start-from-truth and the frozen-tree exact
  gate pass. The same move is dec-B404's candidate for logistic, where
  this survey measured k's 99th percentile from about 23,000 sweeps to
  66.
- **WALNUTS on all leaves and log k** with the latents integrated out
  mixes best per sweep (worst 30 for k, 68 for the fitted values) but a
  sweep costs 27 times as much at the census size and 174 times at 5000
  rows, so per second it loses to the collapsed move everywhere, about
  ties what is built on the slowest census dataset and loses to it
  badly at 5000 rows.
- **Blocking** (k with every leaf integrated out, then all leaves
  jointly) and **the joint leaf draw** condition on the latents, so
  they keep the 1/sqrt(2n) limit: on the slowest datasets the median
  falls by less than half and the worst case gets no better, at 22
  times the cost of a sweep at the census size and 67 to 71 at 5000
  rows. **Calibrated data augmentation**, which needs such a block to be
  exact, gained three to five times on the two slowest datasets once
  tuned (the collapsed move nine to fifteen there), at the same cost.
- A location move and, after the collapsed move, a collapsed location
  move or an elliptical slice on the leaves change little.

**Recommendation.** Before 1.0-0, build the collapsed rescaling move
once, for probit and logistic together, in its lean form (reading the
forest's cached fits and scaling them in place, as the built step
does), taken every fourth sweep beside the built step; for logistic
this is dec-B404's move, so it rides on that study keeping a drawn k.
Do not build blocking, the joint leaf draw, calibrated data
augmentation or the location moves. Do not vendor WALNUTS for this.

## The question

dec-B371's rescaling step (built) took the census's 99th percentile of
k's autocorrelation time from about 23,000 sweeps to about 800, but the
slowest tenth of the small, near-separated datasets still needs several
hundred sweeps per independent draw of k and of the fitted values.
dec-B415 asked for a survey of samplers designed for this kind of
problem - hierarchical-scale funnels and slow data augmentation for
binary responses near separation - prototyped and measured against what
is built, for how much they cut the tail, at what cost, and whether
they are exact.

## Why a tail remains

A probit forest's fit has an overall scale that k, the leaves and the
latent variables share: given the leaves, k is pinned within a few
percent (the funnel); given the latents, the fit's scale is pinned; and
given the fit, the latents sit within about one unit of it. The built
step moves along that shared direction, multiplying the latents and the
leaves by one factor and dividing k by it. Because the latents move
with it, the factor's conditional is set by the latents' residual sum
of squares about the fit, which has n degrees of freedom: the step in
log k is about 1/sqrt(2n) per sweep at every scale of the fit. That is
what reached separation, and it is also the limit: where the posterior
of log k is wide the chain walks it in small steps.

On the census's ten slowest datasets the posterior standard deviation
of log k is 0.70 to 0.94; at n = 150 a step of 0.058 needs some
(0.8 / 0.058)^2, about 200, steps to move one standard deviation, and
the measured autocorrelation times are 400 to 960 sweeps. On the fast
datasets the posterior is 0.26 to 0.40 wide and the times 30 to 100.
At larger n the posterior narrows with the step, so the limit matters
less (below).

Two consequences sort the candidates before any measurement:

- anything that updates k or the leaves given the latents (blocking,
  the joint leaf draw, calibrated augmentation's Gaussian block) can
  remove the funnel and the coupling between trees but keeps the
  latents' grip on the scale, so it cannot cut this tail;
- anything that integrates the latents out of the scale update (the
  collapsed move, WALNUTS or an elliptical slice on the leaves against
  the probit likelihood) takes steps the size of the posterior's own
  width.

## The candidates

Each runs once a sweep ahead of the trees, after the built step, given
the current tree structures.

1. **Blocking.** k drawn from its conditional with every occupied leaf
   integrated out given the latents (a Gaussian marginal whose matrix is
   the total leaf count square, decomposed once a sweep so k's
   conditional costs little to evaluate), then all leaves drawn jointly.
   A partially collapsed Gibbs step (van Dyk and Park 2008).
2. **The joint leaf draw.** All leaves jointly given the latents and k,
   removing backfitting's correlation between trees' leaves.
3. **Calibrated data augmentation** (Duan, Johndrow and Dunson 2018).
   Latents drawn with an inflated variance per row (largest for rows
   the fit already separates), the leaves and k drawn through the block
   of 1 under those latents, the result accepted on the ratio of the
   true to the calibrated likelihood. For that acceptance to be exact
   the update given the latents must be reversible, which a tree-by-tree
   sweep is not; the block supplies it. The variances are set from the
   current fit during burn-in and then fixed.
4. **The collapsed rescaling move.** One factor r: every leaf times r,
   k divided by r, r drawn by one slice step from its conditional with
   the latents integrated out (the probit likelihood at the scaled fit,
   k's prior at k / r), then the latents redrawn given the new fit. The
   form dec-B371's comparison measured and recorded as a door; measured
   here every sweep and every fourth sweep.
5. **WALNUTS** (Bou-Rabee, Carpenter, Kleppe and Liu 2026) on the
   standardized leaves k mu / c and log k, the latents integrated out,
   then the latents redrawn. The standardized form removes the funnel.
   WALNUTS adapts its step within each trajectory, which suits a target
   whose dimension (the leaf count) changes every sweep and so cannot
   carry a tuned mass matrix from sweep to sweep: unit mass, a step of
   0.3, up to six doublings and six halvings.
6. **A location move** (the location half of Zens, Fruhwirth-Schnatter
   and Wagner's boosted samplers for imbalanced binary data): the
   latents and the fit shifted together by one amount, spread evenly
   over every tree's leaves, drawn exactly given the latents.
7. Two additions to 4: **a collapsed location move** (the shift drawn
   against the probit likelihood) and **an elliptical slice** (Murray,
   Adams and MacKay 2010) on all leaves against the probit likelihood.
8. **4 then 1**, to see whether blocking adds to the collapsed move.

## Results

The census: 100 datasets from the probit-k arm's generator (150 rows,
3 predictors, 50 trees, k from chi(1.5, 2)), each drawn by its own
seeded sampler so that every candidate sees the same data and the same
start. With the built step alone these datasets give a median of 73
and a 99th percentile of 879, against the published census's 71 and
about 800. tau is the integrated autocorrelation time in sweeps, of
log k and of the average fitted value.

    all 100 datasets                tau(k): median  90th  99th  worst   tau(fit): median  90th  99th  worst
    built step                                  73   458   879    959                26   550  1003   1050
    + collapsed move (4)                        24    91   124    129                11    31    64     94
    + collapsed move every 4th sweep            35   116   160    166                12    46   112    114
    + 4 and collapsed location (7)              17    58   135    137                11    21    68     77
    + 4 and elliptical slice (7, 60 sets)       20    57    96    120                11    27    52     55
    + location move (6)                         67   418   707    935                22   405   779    869
    + WALNUTS (5), 20 sets *                    13    22    28     30                13    30    64     68

    * the ten slowest datasets under the built step and ten at random

The ten slowest datasets, every candidate (the costly ones on chains of
60,000 recorded sweeps, the rest 200,000):

    candidate                         tau(k): median  worst   tau(fit): median  worst
    built step                                   702    959              642   1050
    blocking (1)                                 393   1679              444   1429
    joint leaf draw (2)                          431   2173              492   1012
    calibrated augmentation (3), 2 sets       188 and 323             197 and 314
    collapsed move (4)                            77    129               12     30
    collapsed move every 4th sweep                79    160               13     45
    4 and collapsed location                      36    137               11     38
    4 and elliptical slice                        35    120               10     20
    location move (6)                            403    935              405    652
    WALNUTS (5)                                   16     30               20     68
    4 then blocking (8)                           26     48               11     21

Cost a sweep relative to the built step, segments of each candidate
interleaved on one chain. "Lean" is the collapsed move as it would be
built; the other prototypes rebuild a table of every row's leaf in
every tree each sweep, about one sweep's work that a build would avoid
for the moves that do not need it.

    candidate                         150 rows, 50 trees   5000 rows, 200 trees
    built step, ms a sweep                    0.047                 2.5
    collapsed move, lean                      1.71                  1.49
    collapsed move, lean, every 4th           1.15                  1.12
    collapsed move, prototype                 2.0                   2.3 to 2.6
    4 and collapsed location                  2.6                   3.0
    4 and elliptical slice                    2.9                   3.4
    location move                             1.3                   1.6
    blocking, joint draw, 4 then 1          22 to 23              67 to 71
    calibrated augmentation                   23                    69
    WALNUTS                                   27                   174

At 5000 rows and 200 trees about 500 leaves are occupied; forming
their matrix and decomposing it is what makes blocking, the joint draw
and calibrated augmentation cost 70 sweeps. WALNUTS spends 50 to 60
gradient evaluations a sweep at the census size and more at 5000 rows,
each a pass over every row and tree.

Time per independent draw of k on the slowest dataset each candidate
met, census size: the built step 45 ms; the collapsed move every sweep 10 ms, every
fourth sweep 9 ms; WALNUTS 38 ms; blocking 1.7 s.

**At larger n** (one dataset per cell, 200 trees, 40,000 sweeps; the
signal scaled to near separation, the fit's sign matching the response
on 90 percent of rows, or left moderate, 76 percent):

    rows   signal     tau(k): built  collapsed    tau(fit): built  collapsed
    1000   separated           201        107                 111         76
    5000   moderate            189         89                  13         14
    5000   separated           274        355                 185        163

Here the posterior of log k is 0.09 to 0.18 wide and the built step's
limit binds less; the collapsed move gains up to two times or nothing,
which is why the recommendation takes it every fourth sweep rather than
every sweep: at 12 percent of a sweep it costs little where it does not
help.

**Exactness.** Every candidate is exact by construction except
calibrated augmentation while its variances adapt (burn-in only).
Measured, the start-from-truth drift of log k over lags 300 to 3000
(a correct sampler sits within |z| 3): the collapsed move z +0.57
(R = 2000) and under logistic +1.23 (R = 2000); with the elliptical
slice +0.65; with the collapsed location -3.09 and +0.99 on two
replicates, -1.34 pooled (R = 4000); the location move +0.14; WALNUTS
-2.88 and +0.23 on two replicates, -1.81 pooled (R = 2000). The
frozen-tree exact gate in quick mode passes for all five, worst |z| 2.9 against its bound of 4.5; its pure
arm's batch-means error, the gate's measure of how fast the small-k
mass is reached, is 0.0031 with the collapsed move, 0.0016 with the
elliptical slice, 0.0024 with the collapsed location and 0.0014 with
WALNUTS, against the built step's 0.0064 (0.0081 with the location
move). The gate's mask arm does not exercise the prototypes, which
decline under a row mask. Taking the collapsed move every fourth sweep
is a fixed schedule of moves that each keep the posterior, so it keeps
it too. Blocking, the joint draw and calibrated augmentation were not
put through either test: their verdict was settled on mixing and cost.

## Why each landed where it did

- **Blocking and the joint draw** gain at the median (on the fast
  datasets k's time falls two to five times: there the funnel and the
  trees' coupling are what is left) but not in the tail, as the
  mechanism says, and they cost 22 to 70 sweeps. Blocking after the
  collapsed move takes the slowest datasets' median from 77 to 26
  sweeps, at 22 times the cost.
- **Calibrated augmentation** was designed for many rows and few
  parameters: there it widens a narrow conditional and the correction
  rejects little. With 100 to 500 leaves moved at once the correction
  bites: at the variances its own rule prescribes the acceptance rate
  was 1 to 5 percent and mixing worse than the built step; scaled down
  ten times it accepted 78 percent and gained three to five times on
  the two slowest datasets, about what blocking gains, at the same
  cost; scaled down a hundred times it accepted 96 percent and gained
  less.
- **The collapsed move** is one dimension along the slow direction with
  the latents out of the way; one slice step moves it about a posterior
  width, so it need not run every sweep. What it leaves (worst 129 to
  160 sweeps) lies in directions it does not move, which WALNUTS and 4
  then blocking reach (16 to 48).
- **WALNUTS** moves every direction at once with the latents out of
  the way, so it mixes best per sweep, but each sweep is 50 or more
  passes over the data.
- **The location moves** help a few datasets where the overall level is
  slow (up to two times) and nothing at the median: the tail here is
  scale, not imbalance.

## Logistic and negative binomial

- **The collapsed move** needs only a likelihood that is cheap to
  evaluate, so it carries over unchanged. On a logistic census (the
  same generator, logistic response) it took k's autocorrelation time
  to a median of 16, 90th percentile 38, 99th 66, worst 96, and the
  fitted values' worst to 28, against logistic without it (69 datasets)
  median 343, 90th 3477, 99th 22,570, worst above 35,000 (lower bounds:
  chains of 200,000 sweeps). Every fourth sweep (logistic has no
  built step beside it): median 30, 90th 61, 99th 115, worst 127, the
  fitted values' worst 59. Start-from-truth under logistic z +1.23. Negative binomial has a cheap likelihood too;
  not measured.
- **WALNUTS and the elliptical slice** need only the likelihood and its
  gradient; they apply to both, at the costs above.
- **The built step** does not apply to logistic as its latents are
  drawn today (Polya-Gamma precisions only). Zens, Fruhwirth-Schnatter
  and Wagner's representation draws a latent utility as well, logistic
  about the fit, with a Polya-Gamma precision PG(2, |utility - fit|);
  given those precisions the utility is Gaussian about the fit, and the
  built step applies with the residual sum of squares weighted by them.
  That would give logistic what probit has now, with probit's tail, at
  the cost of a second Polya-Gamma variate a row. The collapsed move
  does better on both families, so this is recorded, not proposed.
- **Blocking, the joint draw, calibrated augmentation** apply through
  the weighted Gaussian form given the Polya-Gamma precisions, with the
  same limit and cost as under probit.

## What other implementations do

BayesTree, the BART package's probit and logistic fits and bartMachine
fix the leaf scale for binary responses, and stochtree samples it only
when asked; dbarts drawing k by default, as 0.9-34's bart2 did, is the
unusual case, so none of the others meet this funnel by default. Their
binary samplers are Albert and Chib's latents (probit) or a latent
scale mixture (logistic) under the tree moves of Chipman, George and
McCulloch; published work on BART's mixing concerns tree structures
(particle Gibbs, grow-from-root, Pratola's proposals), not the leaf
scale. The general literature's answers to a scale funnel - interweaving
and non-centring (Yu and Meng 2011; Papaspiliopoulos, Roberts and Skold
2007), parameter expansion (Liu and Wu 1999), partial collapsing (van
Dyk and Park 2008) and gradient methods built for funnels (WALNUTS) -
are all represented above or in dec-B371's comparison; and its account
of slow augmentation in binary regression (Johndrow, Smith, Pillai and
Dunson 2019: steps that shrink faster than the posterior) is the
mechanism above, whose remedies are calibration (3), boosting (6 and
the built step) and leaving augmentation for the scale update (4, 5).

## Recommendation

- **Before 1.0-0:** the collapsed rescaling move, built once for probit
  and logistic, in its lean form, every fourth sweep, beside the built
  step under probit. For logistic this is dec-B404's move, conditional
  on its study keeping a drawn k; for probit it is the same code and
  cuts the remaining tail about six times for 12 to 15 percent of a
  sweep. Its exact gate is the built step's frozen-tree gate run with
  both moves, with a logistic arm, and the start-from-truth test under
  both families. How often to take it is the plan's to fix with these
  numbers (every sweep costs 50 to 70 percent for a further 0 to 25
  percent on the tail).
- **After the merge, if wanted:** the elliptical slice after the
  collapsed move (on the slowest datasets k's median 77 to 35, the
  fitted values' worst 30 to 20, little elsewhere, for about 45 percent
  more time than the collapsed move every sweep), only if a smaller
  SBC thin or a further cut is wanted.
- **Not at all:** blocking, the joint leaf draw, calibrated data
  augmentation, the location moves, and WALNUTS for this purpose. On
  vendoring WALNUTS: it works here (exact as far as measured, best per
  sweep, robust to a leaf count that changes every sweep), but it brings
  Eigen, and a sweep with it costs 27 to 174 sweeps without it; the
  leaves have no reason of their own to be updated jointly that would
  pay for that.

## Unsettled

- Negative binomial: the collapsed move applies but was not measured.
- Larger n rests on one dataset per cell; whether every fourth sweep is
  the right frequency at large n (where the move gains up to two times
  or nothing) is for the plan's measurement.
- Costs were measured on a loaded laptop; the bench box was busy
  throughout. Ratios of interleaved segments should carry over; the
  plan should confirm the every-fourth-sweep cost on the box.
- The shipped census script draws each dataset from a sampler the
  previous dataset's chain has advanced, so its per-dataset results
  cannot be compared across builds that draw differently; the survey's
  driver gives each dataset its own sampler. Worth carrying into the
  shipped script when the move is planned.

## Evidence

Runs on an arm64 laptop shared with other work (load 20 to 40 on 10
cores), two at a time; costs are medians of interleaved segments on one
chain (seven rounds for the lean move, three to five for the rest), so
load widens their noise but does not shift one candidate against
another.

**Prototype.** Kept outside the tree as scratch/survey/prototype.patch,
against the bartcore engine of 2026-10-09: the candidates in one source
file, selected by the environment variable DBARTS_SURVEY in a hook after
the built step in Chain::runSweeps (the lean move and its every-k option
inline there); WALNUTS from stan4bart's vendored walnutpie headers
(MIT), compiled against RcppEigen. Drivers and results in
scratch/survey: census.R, truth.R and truth2.R (start-from-truth),
timing.R, largen.R; out/.

**Census.** census.R follows benchmarks/R/probit-k-mixing.R's census
except that each dataset gets its own sampler seeded by its index: the
shipped script reuses one sampler whose generator the previous chain
advances, so a candidate that draws differently would see different
datasets. Sweeps of burn-in and recorded at thin 10, as each run's saved
output records them: on the full census the built step, the collapsed
move every sweep, the location move and the elliptical slice, 30,000 and
200,000; the collapsed location and the every-fourth-sweep run, 2,000
and 200,000; on the ten slowest datasets the every-fourth-sweep run and
the elliptical slice, 2,000 and 200,000; WALNUTS, blocking, the joint
draw, 4 then 1 and calibrated augmentation, 2,000 and 60,000; every
logistic census, 2,000 and 200,000. The ten slowest: datasets 5, 20, 21, 40, 52, 56, 73, 85, 87,
90; the subset of 20 adds 9, 16, 25, 34, 44, 46, 51, 72, 74, 93.
Calibrated augmentation: datasets 20 and 5, adaptation over the first
10,000 sweeps.

**Exact gate.** benchmarks/R/probit-k-scale-exact.R quick, the
prototype library first on the library path.

**References.**
- Albert, J. H. and Chib, S. (1993). Bayesian analysis of binary and
  polychotomous response data. JASA 88, 669-679.
- Bou-Rabee, N., Carpenter, B., Kleppe, T. S. and Liu, S. (2026). The
  within-orbit adaptive leapfrog no-U-turn sampler. JMLR 27(113).
  https://arxiv.org/abs/2506.18746
- Duan, L. L., Johndrow, J. E. and Dunson, D. B. (2018). Scaling up data
  augmentation MCMC via calibration. JMLR 19(64), 1-34.
  https://www.jmlr.org/papers/v19/17-573.html
- Johndrow, J. E., Smith, A., Pillai, N. and Dunson, D. B. (2019). MCMC
  for imbalanced categorical data. JASA 114, 1394-1403.
  https://arxiv.org/abs/1605.05798
- Liu, J. S. and Sabatti, C. (2000). Generalised Gibbs sampler and
  multigrid Monte Carlo for Bayesian computation. Biometrika 87, 353-369.
- Liu, J. S. and Wu, Y. N. (1999). Parameter expansion for data
  augmentation. JASA 94, 1264-1274.
- Murray, I., Adams, R. P. and MacKay, D. J. C. (2010). Elliptical slice
  sampling. AISTATS, 541-548.
- Neal, R. M. (2003). Slice sampling. Annals of Statistics 31, 705-767.
- Papaspiliopoulos, O., Roberts, G. O. and Skold, M. (2007). A general
  framework for the parametrization of hierarchical models. Statistical
  Science 22, 59-73.
- Polson, N. G., Scott, J. G. and Windle, J. (2013). Bayesian inference
  for logistic models using Polya-Gamma latent variables. JASA 108,
  1339-1349.
- van Dyk, D. A. and Park, T. (2008). Partially collapsed Gibbs
  samplers. JASA 103, 790-796.
- Yu, Y. and Meng, X.-L. (2011). To center or not to center: that is not
  the question. JCGS 20, 531-570.
- Zens, G., Fruhwirth-Schnatter, S. and Wagner, H. (2024). Ultimate
  Polya Gamma samplers: efficient MCMC for possibly imbalanced binary
  and categorical data. JASA 119, 2548-2559.
  https://arxiv.org/abs/2011.06898
- stochtree, bart() documentation (sample_sigma2_leaf, default FALSE).
  https://rdrr.io/pkg/stochtree/man/bart.html
- walnutpie (MIT): https://github.com/flatironinstitute/walnutpie
