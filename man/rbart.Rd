\name{rbart}
\alias{rbart_vi}
\alias{plot.rbart}
\alias{fitted.rbart}
\alias{extract.rbart}
\alias{predict.rbart}
\alias{residuals.rbart}
\title{Bayesian Additive Regression Trees with Random Effects}
\description{
  Fits a varying intercept/random effect BART model.
}
\usage{
rbart_vi(
    formula, data, test, subset, weights, offset, offset.test = offset,
    group.by, group.by.test, prior = cauchy,
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
    keepTrees = TRUE, keepCall = TRUE,
    seed = NA_integer_,
    keepSampler = keepTrees,
    \dots)

\S3method{plot}{rbart}(
    x, plquants = c(0.05, 0.95), cols = c('blue', 'black'), \dots)

\S3method{fitted}{rbart}(
    object,
    type = c("ev", "ppd", "bart", "ranef"),
    sample = c("train", "test"),
    \dots)

\S3method{extract}{rbart}(
    object,
    type = c("ev", "ppd", "bart", "ranef", "trees"),
    sample = c("train", "test"),
    combineChains = TRUE,
    \dots)

\S3method{predict}{rbart}(
    object, newdata, group.by, offset,
    type = c("ev", "ppd", "bart", "ranef"),
    combineChains = TRUE,
    \dots)

\S3method{residuals}{rbart}(object, \dots)
}
\arguments{
  \item{group.by}{
    Grouping factor. Can be an integer vector/factor, or a reference to such in \code{data}.
  }
  \item{group.by.test}{
    Grouping factor for test data, of the same type as \code{group.by}. Can be missing.
  }
  \item{prior}{
    A function or symbolic reference to built-in priors. Determines the prior over the standard deviation of the random effects. Supplied functions take two arguments, \code{x} - the standard deviation, and \code{rel.scale} - the standard deviation of the response variable before random effects are fit. Built in priors are \code{cauchy} with a scale of 2.5 times the relative scale and \code{gamma} with a shape of 2.5 and scale of 2.5 times the relative scale.
  }
  \item{n.thin}{
    The number of tree jumps taken for every stored sample, but also the number of samples from the posterior of the standard deviation of the random effects before one is kept.
  }
  \item{formula, data, test, subset, weights, offset, offset.test, sigest, sigdf, sigquant, k, power, base, n.trees, n.samples, n.burn, n.chains, n.threads, combineChains, n.cuts, useQuantiles, keepTrainingFits, printEvery, printCutoffs, verbose, keepTrees, keepCall, seed, keepSampler, \dots}{
    Same as in \code{\link{bart2}}.
  }
  \item{object}{
    A fitted \code{rbart} model.
  }
  \item{newdata}{
    Same as \code{test}, but named to match \code{\link{predict}} generic.
  }
  \item{type}{
    One of \code{"ev"}, \code{"ppd"}, \code{"bart"}, \code{"ranef"}, or \code{"trees"} for the posterior of the expected value, posterior predictive distribution, non-parametric/BART component, random effect, or saved trees respectively. The expected value is the sum of the BART component and the random effects, while the posterior predictive distribution is a response sampled with that mean. To synergize with \code{\link{predict.glm}}, \code{"response"} can be used as a synonym for \code{"value"} and \code{"link"} can be used as a synonym for \code{"bart"}. For additional details on tree extraction, see the corresponding subsection in \code{\link{bart}}.
  }
  \item{sample}{
    One of \code{"train"} or \code{"test"}, referring to the training or tests samples respectively.
  }
  \item{x, plquants, cols}{
    Same as in \code{\link{plot.bart}}.
  }
}
\details{
  Fits a BART model with additive random intercepts, one for each factor level of \code{group.by}. For continuous responses:

\itemize{
  \item \eqn{y_i \sim N(f(x_i) + \alpha_{g[i]}, \sigma^2)}{y_i ~ N(f(x_i) + \alpha_{g[i]}, \sigma^2)}
  \item \eqn{\alpha_j \sim N(0, \tau^2)}{\alpha_j ~ N(0, \tau^2)}.
}

For binary outcomes the response model is changed to \eqn{P(Y_i = 1) = \Phi(f(x_i) + \alpha_{g[i]})}. \eqn{i} indexes observations, \eqn{g[i]} is the group index of observation \eqn{i}, \eqn{f(x)} and \eqn{\sigma_y} come from a BART model, and \eqn{\alpha_j} are the independent and identically distributed random intercepts. Draws from the posterior of \eqn{tau} are made using a slice sampler, with a width dynamically determined by assessing the curvature of the posterior distribution at its mode.

\subsection{Out Of Sample Groups}{
  Predicting random effects for groups not in the training sample is supported by sampling from their posterior predictive distribution, that is a draw is taken from \eqn{p(\alpha \mid y) = \int p(\alpha \mid \tau)p(\tau \mid y)d\alpha}. For out-of-sample groups in the test data, these random effect draws can be kept with the saved object. For those supplied to \code{predict}, they cannot and may change for subsequent calls.
}

\subsection{Generics}{
  See the generics section of \code{\link{bart}}.
}
}
\value{
  An object of class \code{rbart}. Contains all of the same elements of an object of class \code{\link{bart}}, as well as the elements:

  \item{ranef}{
    Samples from the posterior of the random effects. A array/matrix of posterior samples. The \eqn{(k, l, j)} value is the \eqn{l}th draw of the posterior of the random effect for group \eqn{j} (i.e. \eqn{\alpha^*_j}{\alpha*_j}) corresponding to chain \eqn{k}. When \code{n.chains} is one or \code{combineChains} is \code{TRUE}, the result is a collapsed down to a matrix.
  }
  \item{ranef.mean}{
    Posterior mean of random effects, derived by taking mean across group index of samples.
  }
  \item{tau}{
    Matrix of posterior samples of \code{tau}, the standard deviation of the random effects. Dimensions are equal to the number of chains times the numbers of samples unless \code{n.chains} is one or \code{combineChains} is \code{TRUE}.
  }
  \item{\code{first.tau}}{
    Burn-in draws of \code{tau}.
  }
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
                     n.samples = 40L, n.burn = 10L, n.thin = 2L,
                     n.chains = 1L,
                     n.trees = 25L, n.threads = 1L)
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{randomeffects}

