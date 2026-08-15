# Validate A Composed Sampler By Simulation-Based Calibration

Simulation-based calibration (Talts et al. 2018) over a one-sweep step
*you* supply, for a model *you* state as a prior draw and a simulation.
For each replication it draws \\\theta_0\\ from the prior, simulates
data at it, runs the step, and ranks each scalar functional's value at
\\\theta_0\\ among its values over the retained draws. Under a kernel
that targets the posterior the prior and simulation imply, those ranks
are uniform; a systematic departure means the composition does not
sample what it claims to - the defect class in which a posterior *mean*
is reported where a draw belongs. It is a diagnostic, not a fix.

## Usage

``` r
dbartsValidateComposition(drawPrior, simulate, init, step, functionals,
                          n.replications = 200L, n.draws = 200L,
                          n.thin = 30L, n.burn = 200L, alpha = 0.05,
                          seed = NULL)
```

## Arguments

- drawPrior:

  A function of no arguments returning one draw \\\theta_0\\ from the
  prior, in whatever shape `simulate` and `init` expect. Nothing is read
  out of it by name.

- simulate:

  A function of \\\theta_0\\ returning one simulated data set from the
  likelihood the composition assumes.

- init:

  A function of \\\theta_0\\ and the simulated data returning the
  chain's initial state. It **must** start the chain *at* \\\theta_0\\:
  the rank is taken of `functionals(init(theta0, data))`, so a state
  built anywhere else ranks the wrong quantity. The state carries
  whatever `step` and `functionals` need, the data included.

- step:

  A function of the state returning the state after **one** sweep of the
  composed kernel - the whole Gibbs/Metropolis-Hastings scan, not one
  block of it.

- functionals:

  A function of the state returning a numeric vector of scalar summaries
  to rank, the same ones in the same order at every state. They may be
  *derived* - a mean fit, a probability, an indicator - and need not be
  parameters; names, if given, label the report.

- n.replications:

  Number of prior draw / simulation / refit replications. The rank
  resolution and the band both come from this.

- n.draws:

  Number of draws retained per replication, \\L\\; ranks run over \\0\\
  to \\L\\.

- n.thin:

  Sweeps per retained draw. SBC assumes the retained draws are
  independent, so this should exceed the composition's autocorrelation
  time; `1` is right only for a step that draws the posterior outright.

- n.burn:

  Sweeps discarded before the first retained draw. The default is not
  `0`: an `init` that is only approximately at \\\theta_0\\ would
  otherwise be flagged for its early draws. An `init` that is exact can
  use `0`.

- alpha:

  Level of the uniformity band, Bonferroni-corrected within the call
  over the number of functionals ranked, so that a whole call passes
  with probability about \\1 - \alpha\\.

- seed:

  Seeds R's random number stream for the whole run, making it
  reproducible; `NULL` uses the stream as it stands. The band's own
  stream is fixed internally and restored on exit either way, so a call
  never leaves that fix behind.

## Details

**The verdict.** The primary statistic is the maximum absolute
difference between the ranks' empirical cdf and the uniform one, against
a simultaneous band simulated from the null (Talts et al. 2018, fig. 1),
which is already corrected for the multiple looks across the rank grid.
A chi-square goodness of fit and a jittered Kolmogorov-Smirnov test are
reported beside it as secondary signals; the verdict does not hinge on
them.

**Ties.** A rank is the number of draws below \\\theta_0\\, which is
uniform only for an atomless law. Every functional is ranked with the
tie-break instead: draws and \\\theta_0\\ alike carry an independent
uniform tag, so an atom - an indicator, a discrete parameter, a
probability that underflows to zero - contributes a uniform increment
rather than piling every replication onto one rank. With no ties it is
the plain count, exactly.

**The blind spot.** The validator compares the composition against the
model `drawPrior` and `simulate` state *between them*. It cannot detect
a host whose prior draw and simulation disagree with each other: a step
that correctly targets the posterior of the wrong model passes, because
that model is the one being validated against. Nor does it localize a
failure to a block of the sweep; it reports which functionals are
miscalibrated, not why.

## Value

An object of class `dbartsCompositionValidation`, a list with components
`ranks` (the `n.replications` by number-of-functionals integer matrix of
ranks, columns named for the functionals), `L`, `n.replications`,
`n.thin`, `n.burn`, `alpha`, `alpha.adjusted` (the Bonferroni-corrected
level actually used), `mean.rank.target` (\\L/2\\), `pass`, and
`verdicts`, a data frame of one row per functional carrying its mean
rank, ecdf-difference statistic, band, chi-square and KS p-values, and
its `"PASS"` or `"FLAG"` verdict. `print` reports the table.

## References

Talts, S., Betancourt, M., Simpson, D., Vehtari, A., and Gelman, A.
(2018) Validating Bayesian inference algorithms with simulation-based
calibration. *arXiv preprint* arXiv:1804.06788.

## See also

[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md),
the sampler a composition drives, and
[`dbartsAugmentation`](https://vdorie.github.io/dbarts/reference/dbartsAugmentation.md),
the augmentation draws a composed step is built from.

## Examples

``` r
# a conjugate normal-normal one-sweep step: theta ~ N(0, 1) and
# y_i | theta ~ N(theta, 1) give theta | y ~ N(v * sum(y), v), v = 1 / (1 + n)
n <- 5L
v <- 1 / (1 + n)

drawPrior <- function() c(theta = rnorm(1L))
simulateData <- function(theta0) rnorm(n, theta0[["theta"]], 1)
initState <- function(theta0, data)
  list(theta = theta0[["theta"]], mean = v * sum(data), sd = sqrt(v))
functionals <- function(state)
  c(theta = state$theta, positive = as.numeric(state$theta > 0))

# the step draws the posterior outright, so its states are independent and its
# init is exact at theta0: this call needs neither thinning nor burn-in
stepExact <- function(state) {
  state$theta <- rnorm(1L, state$mean, state$sd)
  state
}
dbartsValidateComposition(drawPrior, simulateData, initState, stepExact,
                          functionals, n.replications = 50L, n.draws = 100L,
                          n.thin = 1L, n.burn = 0L, seed = 8L)
#> 
#> Call:
#> dbartsValidateComposition(drawPrior = drawPrior, simulate = simulateData, 
#>     init = initState, step = stepExact, functionals = functionals, 
#>     n.replications = 50L, n.draws = 100L, n.thin = 1L, n.burn = 0L, 
#>     seed = 8L)
#> 
#> 50 replications, 100 draws each (thin 1, burn 0)
#> band alpha 0.05 over 2 functional(s) is 0.025; uniform mean rank 50.0
#> 
#>  functional mean.rank ecdf.diff  band chisq.p  ks.p verdict
#>       theta      49.9    0.0871 0.197   0.828 0.805    PASS
#>    positive      49.4    0.0840 0.197   0.419 0.830    PASS

# the defect class this exists to catch: the posterior MEAN reported where a
# draw belongs. 'theta' flags; 'positive' is coarse enough to survive it, power
# being per functional.
stepMean <- function(state) {
  state$theta <- state$mean
  state
}
dbartsValidateComposition(drawPrior, simulateData, initState, stepMean,
                          functionals, n.replications = 50L, n.draws = 100L,
                          n.thin = 1L, n.burn = 0L, seed = 8L)
#> 
#> Call:
#> dbartsValidateComposition(drawPrior = drawPrior, simulate = simulateData, 
#>     init = initState, step = stepMean, functionals = functionals, 
#>     n.replications = 50L, n.draws = 100L, n.thin = 1L, n.burn = 0L, 
#>     seed = 8L)
#> 
#> 50 replications, 100 draws each (thin 1, burn 0)
#> band alpha 0.05 over 2 functional(s) is 0.025; uniform mean rank 50.0
#> 
#>  functional mean.rank ecdf.diff  band  chisq.p     ks.p verdict
#>       theta        46     0.530 0.197 3.54e-84 1.16e-13    FLAG
#>    positive        48     0.106 0.197 3.26e-01 5.27e-01    PASS
#> 
#> FLAG: the flagged functionals' ranks are not uniform, so the composed 
#> kernel does not target the posterior 'drawPrior' and 'simulate' imply. 
```
