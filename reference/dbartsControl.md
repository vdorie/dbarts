# Discrete Bayesian Additive Regression Trees Sampler Control

Convenience function to create a control object for use with a
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) sampler.

Every integer-valued argument here (`n.samples`, `n.cuts`, `n.burn`,
`n.trees`, `n.chains`, `n.threads`, `n.thin`, `printEvery`,
`printCutoffs`, `categoricalExhaustiveCap`, `testFitParallelCutoff`,
`predictParallelCutoff`, `seed`) refuses a fractional double, naming the
argument, rather than silently truncating it - the same whole-number
rule every other count formal in the package follows, since it too is
coerced through this construction.

## Usage

``` r
dbartsControl(
    verbose = FALSE, keepTrainingFits = TRUE, keepFits = TRUE,
    useQuantiles = FALSE,
    levelGibbs = NA,
    keepTrees = FALSE, storage = c("double", "single"),
    n.samples = NA_integer_,
    n.cuts = 100L, n.burn = 200L, n.trees = 75L, n.chains = 4L,
    n.threads = min(dbarts::guessNumCores(), n.chains), n.thin = 1L, printEvery = 100L,
    printCutoffs = 0L,
    categoricalExhaustiveCap = 10L, testFitParallelCutoff = 65536L,
    predictParallelCutoff = 50000L, sparseDensityThreshold = 0.2,
    proposal.probs = c(
        birth_death = 0.6, swap = 0, change = 0.4, perturb = 0,
        rule_gibbs = 0, birth = 0.5),
    seed = NA_integer_, updateState = TRUE, ...)
```

## Arguments

- ...:

  Not used for new code: the channel that lets a retired argument
  spelling reach a message naming its successor, instead of R's own
  “unused argument” error. Any other name is refused. Removed in dbarts
  1.1-0.

- verbose:

  Logical controlling sampler output to console.

- keepTrainingFits:

  Logical controlling whether or not training fits are returned when the
  sampler runs. These are always computed as part of the fitting
  procedure, so disabling will not substantially impact running time.

- keepFits:

  Logical controlling whether every per-observation channel - training
  fits, test fits, the heteroscedastic variance surface, and per-forest
  fits - is returned when the sampler runs; `keepTrainingFits` is the
  narrower, longstanding switch over the training channel alone,
  unaffected by this one. `FALSE` allocates a per-chain one-draw SCRATCH
  buffer for each opted-out channel instead of an `n.samples`-wide
  array, which a per-draw `callback` (see
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md), and
  [`run`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md))
  can read as each draw is produced without the run ever materializing
  the full array; variable counts and the scalar channels (`sigma`, `k`,
  ...) are kilobytes and are always kept, and `keepTrees` is the
  separate, unaffected recompute-from-saved-trees path. `FALSE` is the
  automatic default of
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) (not
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)) when
  `callback` is supplied and `keepFits` is not itself named. **Three
  warnings apply to `callback`**: a callback that crashes takes down the
  whole R session with no condition to catch; an interrupt cannot land
  while a call is running, so a callback that blocks hangs the session;
  and `keepFits = FALSE` means the sampler's `run()` returns those
  channels as `NULL` (present in the returned list, holding nothing)
  rather than omitting them -
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s packaged
  fit instead OMITS them outright, and its
  `plot`/`extract`/`fitted`/`residuals`/`predict` name `keepFits` in the
  resulting error rather than failing on a bare `NULL`. On a
  `family = "hazard"`/`"hazard.probit"`/`"hazard.logistic"` fit, a
  per-draw `callback`'s draw struct counts person-period-EXPANDED rows,
  not subjects - the mapping back to (subject, period) is done in R and
  is not available to the callback.

- useQuantiles:

  Logical to determine if the empirical quantiles of the columns of
  predictors should be used to determine the tree decision rules. If
  `FALSE`, the rules are spaced uniformly throughout the range of
  covariate values.

- levelGibbs:

  A logical adding one extra Gibbs step per sampler iteration, taken
  before the trees are updated: a constant is added to every occupied
  leaf of each tree, with the constants summing to zero over the trees.
  The sum-of-trees function is therefore unchanged - the fits, the
  residual variance and any latent variables see exactly the state they
  would have - while the individual leaf values move, which can improve
  mixing where the ensemble's overall level is split among many trees.
  The shift is drawn from its exact conditional distribution, so the
  posterior being sampled is the same either way. Three values: `TRUE`
  takes the step every iteration, `FALSE` never takes it, and `NA` (the
  default) takes it for a forest exactly when that forest's tree
  structures are frozen - every structural element of `proposal.probs`
  zero - where the leaf values are the only thing left moving. The
  decision is made afresh each iteration, so a mixture frozen with
  `$setModel` part way through a run switches the step on from there.
  Where structure is being proposed the default takes no step and draws
  exactly the values it drew in previous versions. Applies to forests
  with the default constant leaves (including the monotone constraint
  and every binary, count, or survival family built on them); linear and
  Gaussian-process leaves and the heteroscedastic variance forest ignore
  it. Fixed when the sampler is created.

- keepTrees:

  A logical that determines whether or not trees are cached as they are
  sampled. In all cases, the current state of the sampler is stored as a
  single set of `n.trees`. When `keepTrees` is `TRUE`, a set of
  `n.trees * n.samples` trees are set aside and populated as the sampler
  runs. If the sampler is stopped and restarted, samples proceed from
  the previously stored tree, looping over if necessary. The store costs
  about 24 bytes per node per tree per draw per chain - 23 MB per chain
  at 200 trees, 500 draws, and the eight nodes a tree averages at large
  \\n\\ - and does not grow with the number of observations.

- storage:

  A character string selecting the precision of the internal running
  residual. The default `"double"` stores it in double precision and
  produces bitwise-identical draws to previous versions. `"single"` opts
  the running residual into single (32-bit) precision, halving the
  memory traffic of the dominant per-sweep gather for a speedup that is
  largest where memory bandwidth binds - very large `n` and multi-chain
  runs (leaf draws and reductions remain in double precision). Currently
  supported only for continuous (gaussian) responses with the default
  constant leaves; other models signal an error. The reduced-precision
  path changes the sampled values slightly and so is opt-in.

- n.samples:

  A non-negative integer giving the default number of samples to return
  each time the sampler is run. Generally specified by
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)
  instead, and can be overridden on a per-use basis whenever the sampler
  is
  [`run`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md).
  This is a per-`run()` RETURN count, unaffected by `n.thin` - unlike
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s (and
  [`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md)'s
  `ndpost`) same-named argument, which is a one-shot sweep budget
  divided by thinning; see `bart`'s `n.samples` item for the full
  boundary. `0` is accepted here (and by
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)) - a
  sampler meant to be driven by a host loop's own `run()` calls rather
  than this one's; `bart` and `xbart` both return posterior draws and
  refuse a thinned-to-zero budget instead.

- n.cuts:

  A positive integer or integer vector giving the number of decision
  rules to be used for each given predictor. If of length less than the
  number of predictors, earlier values are recycled. If for any
  predictor more values are specified than are coherent, fewer may be
  used. It does not reach a factor predictor of either kind: an
  unordered factor splits on level subsets and an ordered one at its
  declared level midpoints, a factor's grid following its level table
  rather than a count. See the ‘Decision Rules’ section of
  [`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md) for
  how the rules themselves are placed.

- n.burn:

  A non-negative integer determining how many samples, if any, are
  thrown away at the beginning of a run of the sampler.

- n.trees:

  A positive integer giving the number of trees used in the sum-of-trees
  formulation. Default 75, dbarts's own historical choice; BayesTree's
  and [`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md)'s
  default is 200.

- n.chains:

  A positive integer detailing the number of independent chains for the
  sampler to use.

- n.threads:

  A positive integer giving a total thread budget, distinct from
  `n.chains`: chains run in parallel across it, but tree sampling itself
  uses at most one thread per chain, so within-chain parallelism is not
  shipped and a budget above `n.chains` buys the sampler's own sweep
  nothing - it was built and measured (best case about 3 percent faster
  at four threads on one chain, and slower than serial at two and at
  eight) and closed rather than shipped. The surplus is not wasted
  outright: it still feeds the test-fit pool, used at or above
  `testFitParallelCutoff` test rows, and
  [`predict`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)'s
  own fan-out. Defaults to
  [`guessNumCores`](https://vdorie.github.io/dbarts/reference/guessNumCores.md)
  capped at `n.chains`; a larger, explicit budget is accepted, and
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) and
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) warn once
  per fit, naming both counts, that the excess goes unused for sampling.

- n.thin:

  A positive integer determining how many iterations the MCMC chain
  should jump on the decision trees alone before recording a sample.
  Serves to “thin” the samples against serial correlation. `n.samples`
  are returned regardless of the value of `n.thin`.

- printEvery:

  If `verbose` is `TRUE`, every `printEvery` potential samples (after
  thinning) will issue a verbal statement. Must be a positive integer.

- printCutoffs:

  A non-negative integer specifying how many of the decision rules for a
  variable are printed in verbose mode.

- categoricalExhaustiveCap:

  A positive integer of at least 2 and at most 30 giving the number of
  categories PRESENT at a node up to which the grow-from-root scan
  ([`growFromRoot`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md),
  and [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s
  `n.grow.sweeps`) scores every balanced partition of them, \\2^{P-1} -
  1\\ candidates. Above it the scan scores the \\P - 1\\ sorted prefixes
  instead. This is the one setting here that changes the draws: the
  proposal family above the cap is a different one. Raising it is
  affordable in a narrow band and never beyond - the candidate set
  doubles per level - and cost is not what it protects at its default,
  where enumerating 511 candidates measures the same as emitting 11
  prefixes, the per-proposal cost being the histogram over node members
  that precedes the enumeration. The Metropolis-Hastings birth move
  draws its categorical rules from the prior and never enumerates, so
  this reaches only the grow-from-root path. Fixed when the sampler is
  created.

- testFitParallelCutoff:

  A positive integer giving the number of test rows at or above which a
  chain fits its test matrix across its share of the thread budget
  instead of on its own thread. Time only: the two paths compute in the
  same order and return bit-for-bit the same fits. The default sits past
  the measured crossover on the machine it was measured on - four
  threads win 1.54 times at 65536 rows and 2.97 times at 262144, while
  the serial path already costs 16.5 msec per iteration at 32768 - so a
  caller with large test sets and threads to spare may lower it.
  Requires `n.threads` above `n.chains`; with no surplus budget no chain
  has a pool to route into. Fixed when the sampler is created.

- predictParallelCutoff:

  A positive integer giving the number of tree traversals (rows times
  trees times draws) at or above which
  [`predict`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  spreads a replay over `n.threads` workers instead of running it
  inline. Time only: each worker owns its own output range and nothing
  is summed across workers, so a replay is identical however it is dealt
  out. The default is calibrated - the crossover measures between 3e4
  and 1e5 traversals, the fan-out reaching 1.5 times at 1e5 and
  saturating near 3.5 times on four threads - and replaces an
  uncalibrated estimate of 1e7 that left mid-sized replays serial. Fixed
  when the sampler is created.

- sparseDensityThreshold:

  A number in \\\[0, 1\]\\ giving the nonzero fraction at or below which
  a column of a sparse (`dgCMatrix`) design is stored as a rank bitmap
  rather than as one code per row. Memory against time: at the default a
  rank-stored column costs about 38 percent of a sweep in decode and
  holds about 3.5 times less than the dense layout. The layout changes
  no proposal - the same splits are drawn either way, and the split
  counts are identical - but it is not bitwise: a rank-stored column
  hands a leaf its rows in a different order, so the leaf's sufficient
  statistic reassociates and fits and `sigma` can differ in their last
  bits, by an amount that grows slowly with the number of samples. Which
  side of that trade is worth taking depends on the design and the
  machine, which is why it is a setting. Read once, when the sampler's
  data are built, and fixed thereafter: a change through `setControl` is
  refused by name rather than silently relayouting nothing.

- proposal.probs:

  Named numeric vector or `NULL`, optionally specifying the proposal
  rules and their probabilities. Elements should be `"birth_death"`,
  `"swap"`, `"change"`, `"perturb"` and `"rule_gibbs"` to control tree
  structure proposals, and `"birth"` to give the relative frequency of
  birth/death in the `"birth_death"` step. The five structural
  probabilities must sum to one. All five structural probabilities zero
  is the frozen mixture: no structural proposal is made, the tree
  structures stand where they are, and only the leaf values, `sigma` and
  the family's latents keep being drawn, which is how a fitted forest is
  re-sampled as a fixed basis. Under `dbartsControl`'s default
  `levelGibbs = NA` a frozen forest additionally takes the
  level-shifting Gibbs step each iteration, the leaf values then being
  the only thing left to move. An unnamed `"perturb"` or `"rule_gibbs"`
  is taken as zero and resolved before the rest; an unnamed `"swap"` is
  taken as zero, and a single remaining unnamed element takes the
  residual, so `c(birth_death = 0.7)` is birth/death 0.7, swap 0, change
  0.3, perturb 0, rule_gibbs 0 and `c(birth_death = 0.5, change = 0.4)`
  is swap 0.1; naming only the zero-default moves `"swap"`, `"perturb"`
  and `"rule_gibbs"` leaves the birth/death-versus-change split
  undetermined and is an error. The default is
  `c(birth_death = 0.6, swap = 0, change = 0.4, perturb = 0, rule_gibbs = 0, birth = 0.5)`.
  A `"swap"` element exchanges a parent's split rule with a child's; it
  defaults to zero because at production forest sizes it measures as a
  no-op, but with `n.trees = 1` it is the only move that rotates a rule
  up the tree, so single-tree fits should set it positive (0.1 was the
  historical default). A `"perturb"` element displaces one node's split
  point by a single cut position while keeping its variable and the
  tree's shape; it defaults to zero, and only ordinal (numeric) columns
  can be perturbed. A `"rule_gibbs"` element replaces one nog node's
  rule - a node whose two children are both leaves - with a draw from
  that rule's own full conditional over the available ordinal variables
  and their admissible cuts, so its acceptance is one; it defaults to
  zero, it acts only where the node's own rule is ordinal, and it is
  inert on an all-categorical design. Unlike the four engine settings
  above, this one is not fixed at creation:
  [`setControl`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  accepts a changed mixture between runs, installing it with the priors
  exactly as at creation, and a refused install rolls the stored control
  back.

- seed:

  Random number generator seed. Every chain runs its own generator; the
  seed drives a dedicated generator that in turn hands each chain its
  own seed, leaving R's stream untouched. Seeded results do not depend
  on the thread count, and a single-chain run with a given seed
  reproduces the first chain of a multi-chain run with the same seed. If
  equal to `NA`, chain generators are seeded from R's stream at
  creation, so [`set.seed`](https://rdrr.io/r/base/Random.html)
  beforehand suffices for reproducibility; sampling itself never
  advances R's stream.

- updateState:

  Logical setting the default behavior for many
  [sampler](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  methods with regards to the immediate updating of the cached state of
  the object. A current, cached state is only useful when
  [saving](https://rdrr.io/r/base/save.html)/[loading](https://rdrr.io/r/base/load.html)
  the sampler.

## Engine limits

Eleven fixed values decide what the engine will represent and where it
switches strategy. Each was measured on one machine (arm64 macOS, 10
cores); the four that a workload can profitably move are settings here,
and the rest are limits of the representation or of a proposal family,
with nothing found that moving would buy. A recommended range is given
where a measurement supports one.

|  |  |  |  |  |
|----|----|----|----|----|
| **Limit** | **Where** | **Default** | **Recommended** | **Settable** |
| Cuts per predictor | quantized predictor code | 65533 | up to the cap | no; `n.cuts` above it is refused |
| Categorical levels | quantized predictor code | 65535 | any real factor | no; a wider factor is refused by name |
| Ordered-factor levels | quantized predictor code | 65534 | any real factor | no; one code goes to the cut grid |
| Exact categorical enumeration | grow-from-root scan | 10 present levels | 8 to 14 | yes, `categoricalExhaustiveCap` |
| Leaf-regression columns | [`linear`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md) leaves | 8 | 1 to 8 | no; a fixed stack scratch is sized for it |
| Perturb window | the perturb proposal | 1 grid position | 1 | no; acceptance falls with width and nothing gains |
| Test-fit parallel cutoff | test-set fitting | 65536 rows | 8192 to 65536 | yes, `testFitParallelCutoff` |
| Predict parallel cutoff | [`predict`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md) | 50000 traversals | 3e4 to 1e5 | yes, `predictParallelCutoff` |
| Sparse density threshold | `dgCMatrix` storage | 0.2 | 0.05 to 0.5 | yes, `sparseDensityThreshold` |
| Gaussian-process leaf size | [`gp`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md) leaves | 256 observations | 256 to 512 | yes, `gp(max.leaf.size = )` |
| Person-period expansion | [`hazard`](https://vdorie.github.io/dbarts/reference/dbartsFamilies.md) | 1e7 rows | host-dependent | yes, `hazard(max.rows = )` |

Only the categorical enumeration cap changes what is sampled; the two
parallel cutoffs and the density threshold trade time against memory and
return bit-for-bit the same draws either side of themselves. Every
default above is the value the engine used before it was settable, so a
fit that names none of them reproduces earlier results exactly.

## Value

An object of class `dbartsControl`.

## See also

[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)

## Examples

``` r
## a small, single-chain, reproducible control for embedding a sampler in a
## larger Gibbs loop: one sample per run(), no burn-in, state cached on demand
control <- dbartsControl(n.chains = 1L, n.threads = 1L,
                         n.burn = 0L, n.samples = 1L,
                         n.trees = 25L, seed = 7L,
                         updateState = FALSE)
control
#> An object of class "dbartsControl"
#> Slot "binary":
#> [1] FALSE
#> 
#> Slot "verbose":
#> [1] FALSE
#> 
#> Slot "keepTrainingFits":
#> [1] TRUE
#> 
#> Slot "keepFits":
#> [1] TRUE
#> 
#> Slot "useQuantiles":
#> [1] FALSE
#> 
#> Slot "levelGibbs":
#> [1] NA
#> 
#> Slot "keepTrees":
#> [1] FALSE
#> 
#> Slot "storage":
#> [1] "double"
#> 
#> Slot "n.samples":
#> [1] 1
#> 
#> Slot "n.cuts":
#> [1] 100
#> 
#> Slot "n.burn":
#> [1] 0
#> 
#> Slot "n.trees":
#> [1] 25
#> 
#> Slot "n.chains":
#> [1] 1
#> 
#> Slot "n.threads":
#> [1] 1
#> 
#> Slot "n.thin":
#> [1] 1
#> 
#> Slot "printEvery":
#> [1] 100
#> 
#> Slot "printCutoffs":
#> [1] 0
#> 
#> Slot "categoricalExhaustiveCap":
#> [1] 10
#> 
#> Slot "testFitParallelCutoff":
#> [1] 65536
#> 
#> Slot "predictParallelCutoff":
#> [1] 50000
#> 
#> Slot "sparseDensityThreshold":
#> [1] 0.2
#> 
#> Slot "proposal.probs":
#> birth_death        swap      change     perturb  rule_gibbs       birth 
#>         0.6         0.0         0.4         0.0         0.0         0.5 
#> 
#> Slot "seed":
#> [1] 7
#> 
#> Slot "updateState":
#> [1] FALSE
#> 
#> Slot "call":
#> `NA`()
#> 

n <- 50L
x <- matrix(runif(n * 2L), n, 2L)
y <- x[, 1L] - x[, 2L] + rnorm(n, 0, 0.2)
sampler <- dbarts(y ~ x, control = control)
samples <- sampler$run()
str(samples$train)
#>  num [1:50, 1] -0.00619 -0.6884 0.34457 0.22007 -0.6884 ...
```
