\name{rbart}
\alias{rbart_vi}
\title{Bayesian Additive Regression Trees with Random Effects}
\description{
Fits a varying intercept/random effect BART model.
}
\usage{
rbart_vi(
  formula, data, test, subset, weight, offset, offset.test = offset,
  group.by, prior = cauchy,
  sigest = NA_real_, sigdf = 3.0, sigquant = 0.90,
  k = 2.0,
  power = 2.0, base = 0.95,
  n.trees = 75L,
  n.samples = 1500L, n.burn = 1500L,
  n.chains = 4L, n.threads = min(guessNumCores(), n.chains), combineChains = FALSE,
  n.cuts = 100L, useQuantiles = FALSE,
  n.thin = 5L, keepTrainingFits = TRUE,
  printEvery = 100L, printCutoffs = 0L,
  verbose = TRUE,
  keepTrees = TRUE, keepCall = TRUE, ...)}
\arguments{
  \item{group.by}{
    Grouping factor. Can be an integer vector/factor, or a reference to such in \code{data}.
  }
  \item{prior}{
    A function or symbolic reference to built-in priors. Determines the prior over the standard deviation of 
    the random effects. Supplied functions take two arguments, \code{x} - the standard deviation, and
    \code{rel.scale} - the standard deviation of the response variable before random effects are fit.
    Built in priors are \code{cauchy} with a scale of 2.5 times the relative scale and \code{gamma} with a
    shape of 2.5 and scale of 2.5 times the relative scale.
  }
  \item{n.thin}{
    The number of tree jumps taken for every stored sample, but also the number of samples from the
    posterior of the standard deviation of the random effects before one is kept.
  }
  \item{formula, data, test, subset, weight, offset, offset.test, sigest, sigdf, sigquant,
        k, power, base, n.trees, n.samples, n.burn, n.chains, n.threads, combineChains,
        n.cuts, useQuantiles, keepTrainingFits, printEvery, printCutoffs, verbose,
        keepTrees, keepCall, \dots}{
    Same as in \code{\link{bart2}.}
  }
}
\details{
  rbart best bart.
}
\value{
  Results.
}
\author{
  Vincent Dorie: \email{vdorie@gmail.com}
}
\seealso{
  \code{\link{bart}}, \code{\link{dbarts}}
}
\examples{
f <- function(x) {
    10 * sin(pi * x[,1] * x[,2]) + 20 * (x[,3] - 0.5)^2 +
      10 * x[,4] + 5 * x[,5]
}

set.seed(99)
sigma <- 1.0
n     <- 100

x  <- matrix(runif(n * 10), n, 10)
Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

n.g <- 10
g <- sample(n.g, length(y), replace = TRUE)
sigma.b <- 1.5
b <- rnorm(n.g, 0, sigma.b)

y <- y + b[g]

df <- as.data.frame(x)
colnames(df) <- paste0("x_", seq_len(ncol(x)))
df$y <- y
df$g <- g

## low numbers to reduce run time
rbartFit <- rbart_vi(y ~ . - g, df, group.by = g,
                     n.samples = 40L, n.burn = 10L, n.thin = 2L, n.chains = 1L,
                     n.trees = 25L, n.threads = 1L)
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{randomeffects}

