# Does 1.0-0 compute what 0.9-34 computed?

agent: measurement
rng: neutral (nothing in the package changes; this is a measurement)
budget: one benchmark script, one baseline, this note

## Goal

The equivalence harness carries a nine-scenario record of the classic
engine, taken just before it was deleted. That is the whole of the
quantified cross-engine evidence, and nine scenarios is a narrow base
for the claim the release rests on: that 1.0-0 targets the same
posterior as 0.9-34 wherever both can fit the same model. This note
widens the base to twenty-six scenarios and says where the two
releases agree, where they do not, and why.

## How the comparison is built

Both releases are installed into private libraries and the same script
runs under each. Every scenario is written in the 0.9-x vocabulary -
the BayesTree-spelled fitting function, the sampler object driven by
its own `run` method, the crossvalidation function - so nothing about
the call has to be translated except the one name 1.0-0 moved, and
nothing about the model is left to a default. Every prior and control
setting whose default moved between the two releases is pinned
explicitly on both sides: the number of trees, the burn-in and sample
counts, the number of chains, k, the tree prior's power and base, the
residual prior's degrees of freedom and quantile, the number of cut
points, and the tree-move mixture. That last is the one a careless
comparison would get wrong: 0.9-x proposed birth/death, swap and change
with probability 0.5, 0.1 and 0.4, and 1.0-0 ships 0.6, 0 and 0.4, so a
comparison that left the mixture alone would be measuring the kernel
change and the engine at the same time. The first table below sets the
0.9-x mixture on both sides and so isolates the engine; the second
leaves 1.0-0 at its own new default and so shows the kernel change on
top of it.

Data are fixed within a scenario and only the Monte Carlo seed varies,
across twenty seeds, with a thousand draws kept after five hundred
discarded. The spread across seeds is therefore pure Monte Carlo
variability, and a two-sample Welch statistic formed from the twenty
seeds on each side should behave like a standard normal wherever the
two releases target the same posterior. Each scenario is summarised by
the posterior mean fitted value at twenty-five training rows and every
test row, the posterior mean and standard deviation of the residual
scale, the mean split count for each predictor, and a nonlinear
functional the linear summaries do not pin down - the posterior mean
fraction of fitted values above a fixed cutoff, in and out of sample.
Since a single diverged seed can inflate a variance enough to hide a
shift inside the Welch statistic, disjointness of the two seed ranges
is checked as well; under exchangeability that happens with probability
about one in seventy billion per summary.

## The scenarios

Continuous Friedman; probit; a weighted fit with weights of a half, one
and two; a fit carrying zero weights; a continuous offset; a binary
offset; k drawn under the chi hyperprior; k fixed away from the
default; five cut points and a thousand; quantile cut points; a factor
predictor supplied as a data frame column and, separately, as the
indicator matrix written out by hand; an ordered factor; twenty trees,
two hundred trees, and a single tree; five thousand observations; four
chains; a test set drawn wider than the training range; a sampler
driven around a Gibbs loop with the response and the offset both
replaced between sweeps; a mid-chain predictor swap; crossvalidation
over a two-by-three grid, run over several replications and again over
one; that same crossvalidation written out by hand, fold by fold,
through the ordinary fitting function; and a fit whose predictors have
deliberately unequal numbers of cut points - three continuous against
three binary - over a response that is pure noise.

The last two are not feature tests. The hand-rolled crossvalidation is
the control for the crossvalidation rows: it uses nothing but the
fitting function both releases share, so if the two releases disagree
about a crossvalidated loss but agree about the hand-rolled one, the
disagreement is in the crossvalidation bookkeeping and not in the
sampler. The unequal-cut-point fit is the probe for the tree-structure
posterior: with no signal in the response, nothing in the data prefers
one predictor over another, so the split counts read that posterior
directly, and unequal cut counts are the condition under which the
change move's repaired acceptance ratio can differ from the old one.
Every other scenario here gives every predictor the same number of cut
points.

## Table 1: the engines, at the same tree-move mixture

Twenty-six scenarios, 1500 summaries. Worst summary named where the
maximum falls.

| scenario | summaries | max abs z | worst summary |
| --- | --- | --- | --- |
| friedman | 64 | 2.63 | fitted, test row 1 |
| probit | 62 | 2.92 | split count, predictor 3 |
| weighted | 64 | 3.86 | fitted, test row 22 |
| zeroweights | 64 | 47.28 | sigma mean |
| binaryoffset | 62 | 2.94 | split count, predictor 5 |
| fixedk | 64 | 2.41 | fitted, test row 18 |
| cutssmall | 64 | 2.47 | fitted, train row 23 |
| cutslarge | 64 | 2.47 | fitted, test row 8 |
| usequants | 64 | 3.01 | split count, predictor 4 |
| factorframe | 63 | 2.02 | fitted, train row 16 |
| factorindicators | 63 | 2.02 | fitted, train row 16 |
| orderedfactor | 64 | 3.07 | sigma mean |
| trees20 | 64 | 2.73 | fitted, test row 8 |
| trees200 | 64 | 3.33 | fitted, train row 5 |
| singletree | 64 | 2.91 | fitted, test row 24 |
| largen | 64 | 2.30 | split count, predictor 8 |
| chains4 | 64 | 2.63 | fitted, train row 18 |
| testset | 64 | 2.30 | fitted, test row 23 |
| offset | 64 | 4.20 | sigma mean |
| chik | 66 | 3.61 | sigma mean |
| gibbsloop | 64 | 2.68 | fitted, test row 8 |
| setpredictor | 64 | 2.84 | fitted, train row 11 |
| xbart | 14 | 48.94 | loss at 20 trees, k = 1 |
| xbart1rep | 14 | 19.08 | mean loss |
| cvbyhand | 8 | 1.89 | best cell of the hand-rolled grid |
| mixedcuts | 60 | 6.44 | split count, predictor 4 |

Set aside the two crossvalidation rows, the zero-weight row and the
unequal-cut-point row, whose diagnoses are below, and 1348 summaries
remain over twenty-two scenarios. The largest statistic among them is
4.20; six exceed three and one exceeds four. Under the null, with
twenty seeds a side, one expects about six and about four tenths of one
respectively - the Welch statistic has about thirty-eight degrees of
freedom here, not infinitely many, and its tails are correspondingly
heavier. The observed counts are the expected counts.

## Table 2: 0.9-34 against 1.0-0 at its new default mixture

The same twenty-six scenarios, with 1.0-0 left at birth/death 0.6,
swap 0, change 0.4 and 0.9-34 at its own 0.5, 0.1, 0.4.

| scenario | summaries | max abs z | worst summary |
| --- | --- | --- | --- |
| friedman | 64 | 2.93 | fitted, train row 11 |
| probit | 62 | 2.21 | fitted, train row 25 |
| weighted | 64 | 2.34 | fitted, test row 3 |
| zeroweights | 64 | 58.71 | sigma mean |
| binaryoffset | 62 | 3.12 | fitted, train row 17 |
| fixedk | 64 | 2.69 | fitted, test row 18 |
| cutssmall | 64 | 2.56 | split count, predictor 7 |
| cutslarge | 64 | 2.28 | fitted, train row 9 |
| usequants | 64 | 4.81 | split count, predictor 4 |
| factorframe | 63 | 2.45 | fitted, test row 23 |
| factorindicators | 63 | 2.45 | fitted, test row 23 |
| orderedfactor | 64 | 3.10 | fitted, test row 3 |
| trees20 | 64 | 2.57 | fitted, test row 8 |
| trees200 | 64 | 2.18 | fitted, test row 7 |
| singletree | 64 | 2.18 | fitted, test row 18 |
| largen | 64 | 2.12 | split count, predictor 3 |
| chains4 | 64 | 3.50 | fitted, test row 10 |
| testset | 64 | 2.81 | split count, predictor 6 |
| offset | 64 | 4.80 | fitted, train row 16 |
| chik | 66 | 3.05 | sigma mean |
| gibbsloop | 64 | 3.86 | fitted, test row 6 |
| setpredictor | 64 | 2.82 | fitted, train row 3 |
| xbart | 14 | 51.21 | loss at 20 trees, k = 1 |
| xbart1rep | 14 | 22.33 | mean loss |
| cvbyhand | 8 | 2.44 | 20 trees, k = 4, hand-rolled |
| mixedcuts | 60 | 7.02 | split count, predictor 2 |

Setting aside the same four rows, the largest statistic over the
remaining 1348 summaries is 4.81, with ten above three and two above
four - again close to the six and four tenths the null predicts, and
not distinguishable from Table 1. The two tables share the 0.9-34
recording, so they are not independent readings of the null. Running
the new mixture against the old one inside 1.0-0, so that the engine is
held fixed and only the kernel moves, gives the same verdict directly:
over the same 1500 summaries the largest statistic is 4.03 and only two
exceed four, and the unequal-cut-point row is clean at 2.81. The
mixture change moves the chain, not the distribution it converges to.

## What agrees

Twenty-two of the twenty-six scenarios agree to within Monte Carlo
error, and they agree across the whole range of things a user can ask
for: a continuous response and a binary one, positive weights, an
offset on either scale, a fixed k and a drawn one, a cut grid of five
points and one of a thousand, quantile cut points, a factor reached
either through a data frame or through its own indicator columns, an
ordered factor, one tree and two hundred, five hundred observations and
five thousand, one chain and four, a test set outside the training
range, and the two mutation patterns the sampler object exists for -
swapping the response and the offset between sweeps of an outer Gibbs
loop, and swapping the predictor matrix mid-chain. The factor scenarios
also agree with each other within each release, which says the data
frame and the hand-written indicator matrix reach the same model in
both.

Two engine changes that were expected to show up did not. The initial
forest is now drawn by rejecting whole trees until no leaf is empty
rather than by drawing from the prior and collapsing empty nodes, so
the law at sweep zero differs from every released version of the
package; five hundred discarded draws are enough that nothing of it
survives into the summaries. And the chi hyperprior's degrees-of-
freedom argument was relabelled, so the same written prior means
something different in the two releases; at the infinite scale the
scenario uses, the posterior for k is set by the data rather than by
the hyperprior, and sweeping the degrees of freedom from 1.25 to 5
moves the posterior mean of k by about 0.04 on either side, against a
seed-to-seed spread of 0.07 to 0.09. Neither is contradicted by this
measurement; both are simply too small to see here, which is a
statement about this measurement's resolution and not a proof that they
are negligible everywhere. The third change that was expected to show
up, the change move's repaired acceptance ratio, does show up, and is
the third difference below.

## What differs, and why

Three things differ, and all three are decided changes rather than
surprises.

The first is the zero-weight fit. Carrying rows at weight zero,
0.9-34 reports a posterior mean residual scale of 0.29 where 1.0-0
reports 0.72, and the two seed ranges do not overlap. The anchor
settles which is right: dropping the zero-weight rows and fitting the
remaining ones with no weights at all - a model both releases fit
identically - gives 0.73 under 0.9-34 and 0.73 under 1.0-0. 1.0-0's
weighted fit reproduces that reference; 0.9-34's does not. The
symptom is visible in the fit as well as in the scale: under 0.9-34
the in-sample root mean squared error on the retained rows is about a
tenth, which is near interpolation, against about four tenths under
1.0-0. What 0.9-34 did with a zero-weight row was count it in the
residual-variance posterior's degrees of freedom while letting it
contribute nothing to that posterior's residual sum: the draw divided a
sum over 425 rows by a chi-squared on 500. 1.0-0 counts only rows with
positive weight, so the two agree. That is a defect rather than a
convention, and 0.9-34's own documentation says so - it warns that
weights of zero will be ignored, which is what 1.0-0 does and what
0.9-34 did not. The mismatch is also the whole of the gap. Give those
75 rows a weight of a ten-billionth instead of zero, so that 1.0-0's
own count of positive-weight rows reaches 500, and 1.0-0 collapses to
the same place: residual scale 0.26 against 0.9-34's 0.29, and the same
near-interpolating fit at a root mean squared error of 0.11. An eight
per cent error in the divisor grows into a gap this large because it
compounds - a residual scale drawn too small makes the trees chase
noise, which makes the residual sum smaller, which makes the next draw
smaller still. A second change rides along and is not separately
visible here: empty or zero-weight leaves are now vetoed on the count
of positive-weight members rather than carrying a finite penalty. With
every weight positive the two counts coincide and nothing moves, which
is what the ordinary weighted scenario shows - it is clean in both
tables. So this is a repair, and the measurement quantifies it: on
this data, 0.9-x understated the residual scale by a factor of about
two and a half.

The second is crossvalidation. Every cell of the grid comes out lower
under 0.9-34 than under 1.0-0 - a mean loss of 1.34 against 1.56 - and
of the fourteen summaries the row carries, seven have disjoint seed
ranges, including five of the six grid cells. The hand-rolled control
decides this one too. Run fold by fold through the ordinary fitting
function, the same grid agrees between the releases on every cell, with
a largest statistic of 1.89 over eight summaries. Its values run from
1.37 at the best cell to 1.76 at the worst, which is within three
hundredths of 1.0-0's crossvalidated losses at four of the six cells
and a long way from 0.9-34's at all six: cell by cell the hand-rolled
loss sits 4 to 32 standard errors above 0.9-34's crossvalidated one,
and the whole grid is 0.26 higher on average. So the sampler is not the
difference. The difference is that 0.9-x carried one chain
across the folds and across the replications of the crossvalidation
loop, warm-starting each fold on a forest fitted with that fold's
observations still in the training set, which flatters the reported
loss; 1.0-0 fits each fold fresh. The effect is largest exactly where
the warm start helps most - twenty trees and little shrinkage, where
0.9-34 reports 1.33 against 1.63 - and vanishes at the most-shrunk,
most-averaged cell, a hundred trees at k = 4, where the two agree at
1.34 and 1.37. Two further changes ride along in these rows and are
not separated from it: the third element of the burn-in specification,
a per-replication warm-up, is refused outright now, and each grid cell
is given its own random stream. Neither could produce a one-sided shift
of this size.

The two cells where the hand-rolled control and 1.0-0's crossvalidation
are further apart than that - 0.08 and 0.14 - are both at twenty trees
with little shrinkage, and the cause is not a leak. 1.0-0 runs a fold's
grid cells in one pass, most-shrunk first, so the k = 1 cell inherits
the k = 4 cell's chain and is given only the shorter second burn-in to
forget it. Run that cell on its own, where it gets a fresh sampler and
the full burn-in, and it reports 1.72
against the hand-rolled control's 1.73 to 1.75 - agreement. Nothing
held out enters either fit: with the response randomly permuted, so
that no honest procedure can do better than the response's own spread,
1.0-0's crossvalidation reports 4.97 and the hand-rolled control 4.96.

The third is the change move. Where the predictors have unequal
numbers of cut points, the two releases put different amounts of split
mass on them: over three continuous predictors and three binary ones,
with the response pure noise, 0.9-34 splits on the binary predictors
0.87 times for every split on a continuous one where 1.0-0 splits 0.91
times, a difference of twelve standard errors, and the per-predictor
split counts carry five statistics above four. A fresh block of twenty
seeds gives the same answer, at a largest statistic of 7.2 against 6.4.
This is the change move's repaired acceptance ratio and nothing else,
which three controls show. Take the change move out of the mixture on
both sides and the two releases agree, 0.916 against 0.909, well inside
Monte Carlo error - so none of the other differences between the two
engines produces this. Raise the change move from four tenths of the
proposal mass to eight and the gap more than doubles, from 0.042 to
0.110, which is what a per-move bias does; across the three mixtures
1.0-0's ratio barely moves, at 0.909, 0.909 and 0.910, while 0.9-34's
falls from 0.916 to 0.868 to 0.800. And moving 1.0-0 from the old
mixture to its own new one leaves this row clean, so it is the
acceptance ratio that moved and not the mixture. The old ratio omitted
a proposal-density term that cancels only when the old and new split
variables offer the same number of valid rules,
which is why every other scenario here - all of which give every
predictor the same number of cut points - shows nothing. The size of
the effect is a property of this design, not a general figure: it grows
with the change move's share of the proposal mass and with how unequal
the cut counts are.

## Wall time

Posterior agreement says nothing about speed, and the speed claim the
release makes needs a measurement of both releases on one machine.
This is that measurement, taken 2026-09-14 on the project's
x86 box: an Intel Core i3-14100, four cores, 24 GiB, no swap, running
x86-64 Linux, with a one-minute load average of 0.50 before the run and
1.39 after. 1.0-0 was built from this branch and 0.9-34 from `git
archive main` of this repository, each installed with
`R CMD INSTALL --preclean` into its own private library, and each fit
ran in a fresh `Rscript` process with `R_LIBS` pinned to one of the
two, so only one dbarts ever loaded.

Both releases fit through the BayesTree-style door in 0.9-x's own
argument vocabulary - `bartBT` under 1.0-0, `bart` under 0.9-34 - with
every default that moved between the releases pinned explicitly, the
same pinning the scenarios above use:

    bartFn(x.train = x, y.train = y,
           sigest = NA_real_, sigdf = 3.0, sigquant = 0.90,
           k = 2.0, power = 2.0, base = 0.95, binaryOffset = 0.0,
           ntree = 200L, ndpost = 1000L, nskip = 500L, keepevery = 1L,
           keeptrainfits = TRUE, usequants = FALSE, numcut = 100L,
           verbose = FALSE, nchain = <1|4>, nthread = <1|4>,
           combinechains = TRUE, keeptrees = FALSE, keepcall = FALSE,
           proposalprobs = c(birth_death = 0.5, swap = 0.1,
                             change = 0.4, birth = 0.5))

The design is Friedman data at ten predictors, continuous except in the
probit cell, with the data fixed per cell so only the fit is timed. No
cell carries a factor predictor, so the indicator expansion both
releases do at this door is inert here. Five repetitions per cell,
alternating which library ran first each repetition, so any drift over
the run falls on both alike. Medians of the five:

| cell | design | median 1.0-0 (s) | median 0.9-34 (s) | ratio |
| --- | --- | --- | --- | --- |
| a | n = 1000, one chain | 0.573 | 0.994 | 1.73 |
| b | n = 10000, one chain | 4.155 | 10.615 | 2.55 |
| c | n = 1000, four chains on four threads | 0.659 | 1.130 | 1.71 |
| d | n = 1000, probit, one chain | 0.646 | 1.038 | 1.61 |

The spread within a cell is under two and a half per cent of its median
on both sides, so the ratios are not close calls. The gap widens with
the number of observations, which is what the engine's per-sweep work
predicts: the run loop's running residual and the fused sufficient-
statistic pass both save more where there is more data to sweep.

One machine, one architecture, one design: these are the numbers a user
would quote, not a scaling law.

### xbart k-fold

The TODO entries "xbart follow-through" and "vectorized draw path" record
a regression measured 2026-09-08, before grid cells and folds were
redistributed across workers (landed 2026-09-09): one-replication k-fold
crossvalidation ran 2.1x slower in wall time than 0.9-34 at four threads
while using 33 percent less CPU. This re-measures it on the same x86 box
as the table above - an Intel Core i3-14100, four cores, 24 GiB, no swap,
x86-64 Linux - with a one-minute load average of 0.18 before the run and
1.12 after. 1.0-0 was built from this branch (tip ea8e8e80) and 0.9-34
from `git archive main`, each installed with `R CMD INSTALL --preclean`
into its own private library, each fit run in a fresh `Rscript` process
with `R_LIBS` pinned to one of the two.

The TODO entries do not state a data design, so this uses Friedman data
at ten predictors, continuous, n = 2000, with every moved default pinned
on both sides the way the table above pins them: 200 trees, k = 2,
power = 2, base = 0.95, and `resid.prior = chisq(3.0, 0.90)` (xbart has
no flat `sigdf`/`sigquant`). `xbart` keeps its name and its `control=`
and `n.threads=` spelling on both releases (only the BayesTree-style door
was renamed), so one call runs unmodified under either:

    xbart(y ~ ., data.frame(y = y, x),
          verbose = FALSE, n.samples = 500L,
          method = "k-fold", n.test = 10L, n.reps = 1L,
          n.burn = c(250L, 150L), loss = "rmse",
          n.threads = <1|4>, n.trees = 200L, k = 2, power = 2.0, base = 0.95,
          drop = TRUE, resid.prior = chisq(3.0, 0.90),
          control = dbartsControl(n.cuts = 100L, n.thin = 1L,
                                   updateState = FALSE
                                   [, proposal.probs = c(birth_death = 0.5,
                                                          swap = 0.1, change = 0.4,
                                                          birth = 0.5)]))

10-fold, one replication, 500 draws after 250 burn-in per fold (the
`n.burn` second element, the burn when moving between parameter settings,
is inert here - there is only one grid cell). The bracketed
`proposal.probs` is passed only under 1.0-0, which is the release that
gained the slot (dec-B91); 0.9-34's `dbartsControl` carries no such field
and its own default already equals this, the 0.9-x mixture, so nothing
needs pinning on that side. Five alternating repetitions per cell, wall
and user+system CPU seconds from `proc.time()` around the fit call alone.
Medians of the five:

| n.threads | median 1.0-0 wall (s) | median 0.9-34 wall (s) | ratio | median 1.0-0 CPU (s) | median 0.9-34 CPU (s) |
| --- | --- | --- | --- | --- | --- |
| 1 | 5.513 | 7.397 | 1.34 | 5.510 | 7.390 |
| 4 | 1.896 | 7.382 | 3.89 | 0.045 | 7.380 |

The regression is closed, and reverses: at four threads 1.0-0 now fits
3.89x faster than 0.9-34 in wall time, not 2.1x slower. The mechanism is
exactly what the redistribution changed. This design has a single grid
cell, so 0.9-34's own parallel unit - a (replication, parameter-setting)
pair - never splits below one, and it shows: its wall time does not move
between one thread and four (7.397s to 7.382s). 1.0-0's unit is now the
(fold, grid-cell) pair, ten of them here even at n.reps = 1, so it drops
from 5.513s to 1.896s. The spread within a cell is under 0.8 percent of
its median on both releases at both thread counts, so neither ratio is a
close call.

The four-thread CPU column does not support a "less CPU" claim the way
the 2026-09-08 measurement made one. 1.0-0's reading there (0.045s of a
1.896s fit) does not reflect real CPU-seconds: an external wall/CPU
wrapper around the whole process and system-wide per-core sampling taken
during a run both corroborate the same low number, ruling out an
R-side accounting bug specifically, but a four-thread, compute-bound
control microbenchmark on the same box correctly reports the CPU its
threads use, ruling out a container-wide accounting failure too - and
the fit's own reported loss is bit-identical between one and four
threads, so the same computation runs rather than being skipped. The
likely cause is coarse timer-tick CPU accounting missing short, bursty
per-thread execution on this virtualized box, not anything the engine
does; read the wall-time ratio, not this CPU column, as the four-thread
result. At one thread, where nothing is parallel, CPU tracks wall time
on both releases and the ratio agrees with the wall-time one.

## What stays unmeasured

This is a comparison of what a user gets back from a fit, at the
settings a user is likely to use. It does not reach the sampler's own
accessors - the sums of squared residuals in particular, whose scaling
convention changed, is never read here - nor the state-saving and
tree-printing surfaces, nor prediction from a saved sampler. It does
not reach anything 0.9-34 cannot fit, which is most of what 1.0-0
added: the Student-t, multinomial, ordinal, negative-binomial, hazard
and hurdle families, the two-forest and variance-forest models,
monotonicity and interaction constraints, missingness handling, sparse
predictors, and the DART tree prior. Those have their own exact-
posterior gates and their own equivalence baselines; nothing in this
note speaks to them. Within the scenarios it does cover, the
resolution is set by twenty seeds: flagging at four standard errors
takes a shift of about a third of a posterior standard deviation on a
fitted value, so "agrees" here means "agrees to the precision twenty
seeds can resolve", not "is identical". Finally, everything here was
run on
one machine and one architecture; a platform difference in the
floating-point library would not show up.

The marginal flags are worth a word, because a reader will ask. In
Table 1 a single summary out of 1348 exceeds four, and in Table 2 two
do. Re-running every scenario that carried a statistic above three, on
a fresh block of twenty seeds, cleared both of them - down to 2.2 and
2.5, and 2.3 and 2.4 - and put two different scenarios above four
instead. That is what a nominal-level statistical gate looks like when
it is working, and it is why the harness reports the counts alongside
the maximum rather than gating on the maximum alone. The counts
themselves should be read loosely: the summaries within a scenario are
far from independent of one another, being posterior means of the same
chains, so the spread of these counts around their expectation is
wider than a count of independent tests would give.

## Reproducing

Install 0.9-34 and this release into private libraries, then record
from each and compare the two recordings:

    R_LIBS=<lib-0.9-34> Rscript benchmarks/R/classic-compare.R record old.rds
    R_LIBS=<lib-1.0-0>  Rscript benchmarks/R/classic-compare.R record new.rds mixture=classic
    R_LIBS=<lib-1.0-0>  Rscript benchmarks/R/classic-compare.R record newdef.rds mixture=default
    Rscript benchmarks/R/classic-compare.R compare old.rds new.rds
    Rscript benchmarks/R/classic-compare.R compare old.rds newdef.rds

The 0.9-34 recording is kept in `benchmarks/baselines` as
`classic-compare-0.9-34.rds`, since the release it came from is fixed;
it has to be taken again only if a scenario is added, and it reproduces
bitwise when it is. The 1.0-0 side is the compare side and is not kept.
The wall-time table above is taken separately, with
`benchmarks/R/classic-timing.R`: one cell per invocation under one of the
two libraries, five repetitions alternating which library goes first, the
median of each (cell, library) pair. The driver loop is in
benchmarks/README.md.

`CLASSIC_COMPARE_SCENARIOS` restricts a run to named scenarios,
`CLASSIC_COMPARE_SEED_OFFSET` shifts the seed block for re-checking a
marginal flag, `CLASSIC_COMPARE_CORES` sets the worker count, `merge`
joins chunked recordings, and `quick` runs a smoke-sized version that
is not comparable to a full one. A full recording takes under two
minutes on seven cores.
