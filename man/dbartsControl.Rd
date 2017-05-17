\name{dbartsControl}
\alias{dbartsControl}
\title{Discrete Bayesian Additive Regression Trees Sampler Control}
\description{
  Convenience function to create a control object for use with a \code{\link{dbarts}} sampler.
}
\usage{
dbartsControl(verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
              n.samples = NA_integer_, n.cuts = 100L,
              n.burn = 200L, n.trees = 75L, n.chains = 4L, n.threads = guessNumCores(),
              n.thin = 1L, printEvery = 100L, printCutoffs = 0L,
              rngKind = "default", rngNormalKind = "default", updateState = TRUE)
}
\arguments{
   \item{verbose}{Logical controlling sampler output to console.}
   \item{keepTrainingFits}{Logical controlling whether or not training fits are returned when
     the sampler runs. These are always computed as part of the fitting procedure, so disabling
     will not substantially impact running time.}
   \item{useQuantiles}{Logical to determine if the empirical quantiles of a columns of predictors
   	 should be used to determine the tree decision rules. If \code{FALSE}, the rules are spaced
   	 uniformly throughout the range of covariate values.}
   \item{n.samples}{A non-negative integer giving the default number of samples to return each time the
   	 sampler is run. Generally specified by \code{\link{dbarts}} instead, and can be overridden
   	 on a per-use basis whenever the sampler is \code{\link[=dbartsSampler-class]{run}}.}
   \item{n.cuts}{A positive integer or integer vector giving the number of decision rules to be used
   	 for each given predictor. If of length less than the number of predictors, earlier values
   	 are recycled. If for any predictor more values are specified than are coherent, fewer may
   	 be used. See details for more information.}
   \item{n.burn}{A non-negative integer determining how many samples, if any, are thrown away at the beginning
   	 of a run of the sampler.}
   \item{n.trees}{A positive integer giving the number of trees used in the sum-of-trees formulation.}
   \item{n.chains}{A positive integer detailing the number of independent chains for the sampler to use.}
   \item{n.threads}{A positive integer controlling how many threads will be used for various internal
   	 calculations, as well as the number of chains. Internal calculations are highly optimized so that
         single-threaded performance tends to be superior unless the number of observations is very large
         (>10k), so that it is often not necessary to have the number of threads exceed the number of chains.}
   \item{n.thin}{A positive integer determining how many iterations the MCMC chain should jump on the
   	 decision trees alone before recording a sample. Serves to \dQuote{thin} the samples against
   	 serial correlation. \code{n.samples} are returned regardless of the value of \code{n.thin}.}
   \item{printEvery}{If \code{verbose} is \code{TRUE}, every \code{printEvery} potential samples
   	 (after thinning) will issue a verbal statement. Must be a positive integer.}
   \item{printCutoffs}{A non-negative integer specifying how many of the decision rules for a variable are
   	 printed in verbose mode.}
   \item{rngKind}{Random number generator kind, as used in \code{\link[=Random]{set.seed}}. For type
         \code{"default"}, the built-in generator will be used if possible. Otherwise, will attempt to
         match the built-in generator's type. Success depends on the number of threads.}
   \item{rngNormalKind}{Random number generator normal kind, as used in \code{\link[=Random]{set.seed}}.
         For type \code{"default"}, the built-in generator will be used if possible. Otherwise, will
         attempt to match the built-in generator's type. Success depends on the number of threads and
         the \code{rngKind}.}
   \item{updateState}{Logical setting the default behavior for many \link[=dbartsSampler-class]{sampler} methods
   	 with regards to the immediate updating of the cached state of the object. A current, cached state
   	 is only useful when \link[=save]{saving}/\link[=load]{loading} the sampler.}
}
\value{
  An object of class \code{dbartControl}.
}
\seealso{
  \code{\link{dbarts}}
}
