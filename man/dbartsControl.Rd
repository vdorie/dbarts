\name{dbartsControl}
\alias{dbartsControl}
\title{Discrete Bayesian Additive Regression Trees Sampler Control}
\description{
  Convenience function to create a control object for use with a \code{\link{dbarts}} sampler.
}
\usage{
dbartsControl(
    verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
    keepTrees = FALSE, storage = c("double", "single"),
    n.samples = NA_integer_,
    n.cuts = 100L, n.burn = 200L, n.trees = 75L, n.chains = 4L,
    n.threads = dbarts::guessNumCores(), n.thin = 1L, printEvery = 100L,
    printCutoffs = 0L, rngSeed = NA_integer_, updateState = TRUE)
}
\arguments{
   \item{verbose}{
     Logical controlling sampler output to console.
   }
   \item{keepTrainingFits}{
     Logical controlling whether or not training fits are returned when the sampler runs. These are always computed as part of the fitting procedure, so disabling will not substantially impact running time.
   }
   \item{useQuantiles}{
     Logical to determine if the empirical quantiles of the columns of predictors should be used to determine the tree decision rules. If \code{FALSE}, the rules are spaced uniformly throughout the range of covariate values.
   }
   \item{keepTrees}{
     A logical that determines whether or not trees are cached as they are sampled. In all cases, the current state of the sampler is stored as a single set of \code{n.trees}. When \code{keepTrees} is \code{TRUE}, a set of \code{n.trees * n.samples} trees are set aside and populated as the sampler runs. If the sampler is stopped and restarted, samples proceed from the previously stored tree, looping over if necessary.
   }
   \item{storage}{
     A character string selecting the precision of the internal running residual. The default \code{"double"} stores it in double precision and produces bitwise-identical draws to previous versions. \code{"single"} opts the running residual into single (32-bit) precision, halving the memory traffic of the dominant per-sweep gather for a speedup that is largest where memory bandwidth binds - very large \code{n} and multi-chain runs (leaf draws and reductions remain in double precision). Currently supported only for continuous (gaussian) responses with the default constant leaves; other models signal an error. The reduced-precision path changes the sampled values slightly and so is opt-in.
   }
   \item{n.samples}{
     A non-negative integer giving the default number of samples to return each time the sampler is run. Generally specified by \code{\link{dbarts}} instead, and can be overridden on a per-use basis whenever the sampler is \code{\link[=dbartsSampler-class]{run}}.
   }
   \item{n.cuts}{
     A positive integer or integer vector giving the number of decision rules to be used for each given predictor. If of length less than the number of predictors, earlier values are recycled. If for any predictor more values are specified than are coherent, fewer may be used. See the \sQuote{Decision Rules} section of \code{\link{bart}} for how the rules themselves are placed.
   }
   \item{n.burn}{
     A non-negative integer determining how many samples, if any, are thrown away at the beginning of a run of the sampler.
   }
   \item{n.trees}{
     A positive integer giving the number of trees used in the sum-of-trees formulation. Default 75, dbarts's own historical choice; BayesTree's and \code{\link{bart}}'s default is 200.
   }
   \item{n.chains}{
     A positive integer detailing the number of independent chains for the sampler to use.
   }
   \item{n.threads}{
     A positive integer controlling how many threads will be used for various internal calculations, as well as the number of chains. Internal calculations are highly optimized so that single-threaded performance tends to be superior unless the number of observations is very large (>10k), so that it is often not necessary to have the number of threads exceed the number of chains.
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
   \item{rngSeed}{
     Random number generator seed. Every chain runs its own generator; the seed drives a dedicated generator that in turn hands each chain its own seed, leaving R's stream untouched. Seeded results do not depend on the thread count, and a single-chain run with a given seed reproduces the first chain of a multi-chain run with the same seed. If equal to \code{NA}, chain generators are seeded from R's stream at creation, so \code{\link[=Random]{set.seed}} beforehand suffices for reproducibility; sampling itself never advances R's stream.
   }
   \item{updateState}{
     Logical setting the default behavior for many \link[=dbartsSampler-class]{sampler} methods with regards to the immediate updating of the cached state of the object. A current, cached state is only useful when \link[=save]{saving}/\link[=load]{loading} the sampler.
   }
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
                         n.trees = 25L, rngSeed = 7L,
                         updateState = FALSE)
control

n <- 50L
x <- matrix(runif(n * 2L), n, 2L)
y <- x[, 1L] - x[, 2L] + rnorm(n, 0, 0.2)
sampler <- dbarts(y ~ x, control = control)
samples <- sampler$run()
str(samples$train)
}
