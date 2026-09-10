\name{dbartsControl}
\alias{dbartsControl}
\title{Discrete Bayesian Additive Regression Trees Sampler Control}
\description{
  Convenience function to create a control object for use with a \code{\link{dbarts}} sampler.

  Every integer-valued argument here (\code{n.samples}, \code{n.cuts}, \code{n.burn}, \code{n.trees}, \code{n.chains}, \code{n.threads}, \code{n.thin}, \code{printEvery}, \code{printCutoffs}, \code{categoricalExhaustiveCap}, \code{testFitParallelCutoff}, \code{predictParallelCutoff}, \code{seed}) refuses a fractional double, naming the argument, rather than silently truncating it - the same whole-number rule every other count formal in the package follows, since it too is coerced through this construction.
}
\usage{
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
    seed = NA_integer_, updateState = TRUE, \dots)
}
\arguments{
  \item{\dots}{
    Not used for new code: the channel that lets a retired argument spelling reach a message naming its successor, instead of R's own \dQuote{unused argument} error. Any other name is refused. Removed in dbarts 1.1-0.
  }
   \item{verbose}{
     Logical controlling sampler output to console.
   }
   \item{keepTrainingFits}{
     Logical controlling whether or not training fits are returned when the sampler runs. These are always computed as part of the fitting procedure, so disabling will not substantially impact running time.
   }
   \item{keepFits}{
     Logical controlling whether every per-observation channel - training fits, test fits, the heteroscedastic variance surface, and per-forest fits - is returned when the sampler runs; \code{keepTrainingFits} is the narrower, longstanding switch over the training channel alone, unaffected by this one. \code{FALSE} allocates a per-chain one-draw SCRATCH buffer for each opted-out channel instead of an \code{n.samples}-wide array, which a per-draw \code{callback} (see \code{\link{bart}}, \code{\link{dbarts}}, and \code{\link[=dbartsSampler-class]{run}}) can read as each draw is produced without the run ever materializing the full array; variable counts and the scalar channels (\code{sigma}, \code{k}, ...) are kilobytes and are always kept, and \code{keepTrees} is the separate, unaffected recompute-from-saved-trees path. \code{FALSE} is the automatic default of \code{\link{bart}} (not \code{\link{dbarts}}) when \code{callback} is supplied and \code{keepFits} is not itself named. \strong{Three warnings apply to \code{callback}}: a callback that crashes takes down the whole R session with no condition to catch; an interrupt cannot land while a call is running, so a callback that blocks hangs the session; and \code{keepFits = FALSE} means the sampler's \code{run()} returns those channels as \code{NULL} (present in the returned list, holding nothing) rather than omitting them - \code{\link{bart}}'s packaged fit instead OMITS them outright, and its \code{plot}/\code{extract}/\code{fitted}/\code{residuals}/\code{predict} name \code{keepFits} in the resulting error rather than failing on a bare \code{NULL}. On a \code{family = "hazard"}/\code{"hazard.probit"}/\code{"hazard.logistic"} fit, a per-draw \code{callback}'s draw struct counts person-period-EXPANDED rows, not subjects - the mapping back to (subject, period) is done in R and is not available to the callback.
   }
   \item{useQuantiles}{
     Logical to determine if the empirical quantiles of the columns of predictors should be used to determine the tree decision rules. If \code{FALSE}, the rules are spaced uniformly throughout the range of covariate values.
   }
   \item{levelGibbs}{
     A logical adding one extra Gibbs step per sampler iteration, taken before the trees are updated: a constant is added to every occupied leaf of each tree, with the constants summing to zero over the trees. The sum-of-trees function is therefore unchanged - the fits, the residual variance and any latent variables see exactly the state they would have - while the individual leaf values move, which can improve mixing where the ensemble's overall level is split among many trees. The shift is drawn from its exact conditional distribution, so the posterior being sampled is the same either way. Three values: \code{TRUE} takes the step every iteration, \code{FALSE} never takes it, and \code{NA} (the default) takes it for a forest exactly when that forest's tree structures are frozen - every structural element of \code{proposal.probs} zero - where the leaf values are the only thing left moving. The decision is made afresh each iteration, so a mixture frozen with \code{$setModel} part way through a run switches the step on from there. Where structure is being proposed the default takes no step and draws exactly the values it drew in previous versions. Applies to forests with the default constant leaves (including the monotone constraint and every binary, count, or survival family built on them); linear and Gaussian-process leaves and the heteroscedastic variance forest ignore it. Fixed when the sampler is created.
   }
   \item{keepTrees}{
     A logical that determines whether or not trees are cached as they are sampled. In all cases, the current state of the sampler is stored as a single set of \code{n.trees}. When \code{keepTrees} is \code{TRUE}, a set of \code{n.trees * n.samples} trees are set aside and populated as the sampler runs. If the sampler is stopped and restarted, samples proceed from the previously stored tree, looping over if necessary. The store costs about 24 bytes per node per tree per draw per chain - 23 MB per chain at 200 trees, 500 draws, and the eight nodes a tree averages at large \eqn{n}{n} - and does not grow with the number of observations.
   }
   \item{storage}{
     A character string selecting the precision of the internal running residual. The default \code{"double"} stores it in double precision and produces bitwise-identical draws to previous versions. \code{"single"} opts the running residual into single (32-bit) precision, halving the memory traffic of the dominant per-sweep gather for a speedup that is largest where memory bandwidth binds - very large \code{n} and multi-chain runs (leaf draws and reductions remain in double precision). Currently supported only for continuous (gaussian) responses with the default constant leaves; other models signal an error. The reduced-precision path changes the sampled values slightly and so is opt-in.
   }
   \item{n.samples}{
     A non-negative integer giving the default number of samples to return each time the sampler is run. Generally specified by \code{\link{dbarts}} instead, and can be overridden on a per-use basis whenever the sampler is \code{\link[=dbartsSampler-class]{run}}. This is a per-\code{run()} RETURN count, unaffected by \code{n.thin} - unlike \code{\link{bart}}'s (and \code{\link{bartBT}}'s \code{ndpost}) same-named argument, which is a one-shot sweep budget divided by thinning; see \code{bart}'s \code{n.samples} item for the full boundary. \code{0} is accepted here (and by \code{\link{dbarts}}) - a sampler meant to be driven by a host loop's own \code{run()} calls rather than this one's; \code{bart} and \code{xbart} both return posterior draws and refuse a thinned-to-zero budget instead.
   }
   \item{n.cuts}{
     A positive integer or integer vector giving the number of decision rules to be used for each given predictor. If of length less than the number of predictors, earlier values are recycled. If for any predictor more values are specified than are coherent, fewer may be used. It does not reach a factor predictor of either kind: an unordered factor splits on level subsets and an ordered one at its declared level midpoints, a factor's grid following its level table rather than a count. See the \sQuote{Decision Rules} section of \code{\link{bartBT}} for how the rules themselves are placed.
   }
   \item{n.burn}{
     A non-negative integer determining how many samples, if any, are thrown away at the beginning of a run of the sampler.
   }
   \item{n.trees}{
     A positive integer giving the number of trees used in the sum-of-trees formulation. Default 75, dbarts's own historical choice; BayesTree's and \code{\link{bartBT}}'s default is 200.
   }
   \item{n.chains}{
     A positive integer detailing the number of independent chains for the sampler to use.
   }
   \item{n.threads}{
     A positive integer giving a total thread budget, distinct from \code{n.chains}: chains run in parallel across it, but tree sampling itself uses at most one thread per chain, so within-chain parallelism is not shipped and a budget above \code{n.chains} buys the sampler's own sweep nothing - it was built and measured (best case about 3 percent faster at four threads on one chain, and slower than serial at two and at eight) and closed rather than shipped. The surplus is not wasted outright: it still feeds the test-fit pool, used at or above \code{testFitParallelCutoff} test rows, and \code{\link[=dbartsSampler-class]{predict}}'s own fan-out. Defaults to \code{\link{guessNumCores}} capped at \code{n.chains}; a larger, explicit budget is accepted, and \code{\link{dbarts}} and \code{\link{bart}} warn once per fit, naming both counts, that the excess goes unused for sampling.
   }
   \item{n.thin}{
     A positive integer determining how many iterations the MCMC chain should jump on the decision trees alone before recording a sample. Serves to \dQuote{thin} the samples against serial correlation. \code{n.samples} are returned regardless of the value of \code{n.thin}.
   }
   \item{printEvery}{
     If \code{verbose} is \code{TRUE}, every \code{printEvery} potential samples (after thinning) will issue a verbal statement. Must be a positive integer.
   }
   \item{printCutoffs}{
     A non-negative integer specifying how many of the decision rules for a variable are printed in verbose mode.
   }
   \item{categoricalExhaustiveCap}{
     A positive integer of at least 2 and at most 30 giving the number of categories PRESENT at a node up to which the grow-from-root scan (\code{\link[=dbartsSampler-class]{growFromRoot}}, and \code{\link{bart}}'s \code{n.grow.sweeps}) scores every balanced partition of them, \eqn{2^{P-1} - 1}{2^(P-1) - 1} candidates. Above it the scan scores the \eqn{P - 1}{P - 1} sorted prefixes instead. This is the one setting here that changes the draws: the proposal family above the cap is a different one. Raising it is affordable in a narrow band and never beyond - the candidate set doubles per level - and cost is not what it protects at its default, where enumerating 511 candidates measures the same as emitting 11 prefixes, the per-proposal cost being the histogram over node members that precedes the enumeration. The Metropolis-Hastings birth move draws its categorical rules from the prior and never enumerates, so this reaches only the grow-from-root path. Fixed when the sampler is created.
   }
   \item{testFitParallelCutoff}{
     A positive integer giving the number of test rows at or above which a chain fits its test matrix across its share of the thread budget instead of on its own thread. Time only: the two paths compute in the same order and return bit-for-bit the same fits. The default sits past the measured crossover on the machine it was measured on - four threads win 1.54 times at 65536 rows and 2.97 times at 262144, while the serial path already costs 16.5 msec per iteration at 32768 - so a caller with large test sets and threads to spare may lower it. Requires \code{n.threads} above \code{n.chains}; with no surplus budget no chain has a pool to route into. Fixed when the sampler is created.
   }
   \item{predictParallelCutoff}{
     A positive integer giving the number of tree traversals (rows times trees times draws) at or above which \code{\link[=dbartsSampler-class]{predict}} spreads a replay over \code{n.threads} workers instead of running it inline. Time only: each worker owns its own output range and nothing is summed across workers, so a replay is identical however it is dealt out. The default is calibrated - the crossover measures between 3e4 and 1e5 traversals, the fan-out reaching 1.5 times at 1e5 and saturating near 3.5 times on four threads - and replaces an uncalibrated estimate of 1e7 that left mid-sized replays serial. Fixed when the sampler is created.
   }
   \item{sparseDensityThreshold}{
     A number in \eqn{[0, 1]}{[0, 1]} giving the nonzero fraction at or below which a column of a sparse (\code{dgCMatrix}) design is stored as a rank bitmap rather than as one code per row. Memory against time, with no effect on the answers: at the default a rank-stored column costs about 38 percent of a sweep in decode and holds about 3.5 times less than the dense layout. Which side of that trade is worth taking depends on the design and the machine, which is why it is a setting. Read once, when the sampler's data are built, and fixed thereafter: a change through \code{setControl} is refused by name rather than silently relayouting nothing.
   }
   \item{seed}{
     Random number generator seed. Every chain runs its own generator; the seed drives a dedicated generator that in turn hands each chain its own seed, leaving R's stream untouched. Seeded results do not depend on the thread count, and a single-chain run with a given seed reproduces the first chain of a multi-chain run with the same seed. If equal to \code{NA}, chain generators are seeded from R's stream at creation, so \code{\link[=Random]{set.seed}} beforehand suffices for reproducibility; sampling itself never advances R's stream.
   }
   \item{updateState}{
     Logical setting the default behavior for many \link[=dbartsSampler-class]{sampler} methods with regards to the immediate updating of the cached state of the object. A current, cached state is only useful when \link[=save]{saving}/\link[=load]{loading} the sampler.
   }
}
\section{Engine limits}{
  Eleven fixed values decide what the engine will represent and where it
  switches strategy. Each was measured on one machine (arm64 macOS, 10 cores);
  the four that a workload can profitably move are settings here, and the rest
  are limits of the representation or of a proposal family, with nothing found
  that moving would buy. A recommended range is given where a measurement
  supports one.

  \tabular{lllll}{
    \strong{Limit} \tab \strong{Where} \tab \strong{Default} \tab \strong{Recommended} \tab \strong{Settable} \cr
    Cuts per predictor \tab quantized predictor code \tab 65533 \tab up to the cap \tab no; \code{n.cuts} above it is refused \cr
    Categorical levels \tab quantized predictor code \tab 65535 \tab any real factor \tab no; a wider factor is refused by name \cr
    Ordered-factor levels \tab quantized predictor code \tab 65534 \tab any real factor \tab no; one code goes to the cut grid \cr
    Exact categorical enumeration \tab grow-from-root scan \tab 10 present levels \tab 8 to 14 \tab yes, \code{categoricalExhaustiveCap} \cr
    Leaf-regression columns \tab \code{\link[=dbartsPriors]{linear}} leaves \tab 8 \tab 1 to 8 \tab no; a fixed stack scratch is sized for it \cr
    Perturb window \tab the perturb proposal \tab 1 grid position \tab 1 \tab no; acceptance falls with width and nothing gains \cr
    Test-fit parallel cutoff \tab test-set fitting \tab 65536 rows \tab 8192 to 65536 \tab yes, \code{testFitParallelCutoff} \cr
    Predict parallel cutoff \tab \code{\link[=dbartsSampler-class]{predict}} \tab 50000 traversals \tab 3e4 to 1e5 \tab yes, \code{predictParallelCutoff} \cr
    Sparse density threshold \tab \code{dgCMatrix} storage \tab 0.2 \tab 0.05 to 0.5 \tab yes, \code{sparseDensityThreshold} \cr
    Gaussian-process leaf size \tab \code{\link[=dbartsPriors]{gp}} leaves \tab 256 observations \tab 256 to 512 \tab yes, \code{gp(max.leaf.size = )} \cr
    Person-period expansion \tab \code{\link[=dbartsFamilies]{hazard}} \tab 1e7 rows \tab host-dependent \tab yes, \code{hazard(max.rows = )} \cr
  }

  Only the categorical enumeration cap changes what is sampled; the two
  parallel cutoffs and the density threshold trade time against memory and
  return bit-for-bit the same draws either side of themselves. Every default
  above is the value the engine used before it was settable, so a fit that
  names none of them reproduces earlier results exactly.
}
\value{
  An object of class \code{dbartsControl}.
}
\seealso{
  \code{\link{dbarts}}
}
\examples{
## a small, single-chain, reproducible control for embedding a sampler in a
## larger Gibbs loop: one sample per run(), no burn-in, state cached on demand
control <- dbartsControl(n.chains = 1L, n.threads = 1L,
                         n.burn = 0L, n.samples = 1L,
                         n.trees = 25L, seed = 7L,
                         updateState = FALSE)
control

n <- 50L
x <- matrix(runif(n * 2L), n, 2L)
y <- x[, 1L] - x[, 2L] + rnorm(n, 0, 0.2)
sampler <- dbarts(y ~ x, control = control)
samples <- sampler$run()
str(samples$train)
}
