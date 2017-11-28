\name{bart}
\alias{bart}
\alias{bart2}
\alias{plot.bart}
\title{Bayesian Additive Regression Trees}
\description{
BART is a Bayesian \dQuote{sum-of-trees} model in which each tree is constrained by a prior to be a weak learner.

\itemize{
  \item For numeric response \eqn{y = f(x) + \epsilon}{y = f(x) + \epsilon},
where \eqn{\epsilon \sim N(0, \sigma^2)}{\epsilon ~ N(0, \sigma^2)}.
  \item For binary response \eqn{y}, \eqn{P(Y = 1 \mid x) = \Phi(f(x))}{P(Y = 1 | x) = \Phi(f(x))}, where \eqn{\Phi}
denotes the standard normal cdf (probit link).
}
}
\usage{
bart(x.train, y.train, x.test = matrix(0.0, 0, 0),
     sigest = NA, sigdf = 3, sigquant = 0.90,
     k = 2.0,
     power = 2.0, base = 0.95,
     binaryOffset = 0.0, weights = NULL,
     ntree = 200,
     ndpost = 1000, nskip = 100,
     printevery = 100, keepevery = 1, keeptrainfits = TRUE,
     usequants = FALSE, numcut = 100, printcutoffs = 0,
     verbose = TRUE, nchain = 1, nthread = 1, keepcall = TRUE)

bart2(formula, data, test, subset, weights, offset, offset.test = offset,
      sigest = NA_real_, sigdf = 3.0, sigquant = 0.90,
      k = 2.0,
      power = 2.0, base = 0.95,
      n.trees = 75L,
      n.samples = 500L, n.burn = 500L,
      n.chains = 4L, n.threads = min(guessNumCores(), n.chains),
      n.cuts = 100L, useQuantiles = FALSE,
      n.thin = 1L, keepTrainingFits = TRUE,
      printEvery = 100L, printCutoffs = 0L,
      verbose = TRUE, keepCall = TRUE, ...)
     
\method{plot}{bart}(x,
     plquants = c(0.05, 0.95), cols = c('blue', 'black'),
     \dots)
}
\arguments{
   \item{x.train}{
     Explanatory variables for training (in sample) data. May be a matrix or a data frame,
     with rows corresponding to observations and columns to variables.
     If a variable is a factor in a data frame, it is replaced with dummies.
     Note that \eqn{q} dummies are created if \eqn{q > 2} and
     one dummy is created if \eqn{q = 2}, where \eqn{q} is the number of levels of the factor.
   }

   \item{y.train}{
     Dependent variable for training (in sample) data. If \code{y.train} is numeric a continous response
     model is fit (normal errors). If \code{y.train} is a binary factor or has only values 0 and 1, then
     a binary response model with a probit link is fit.
   }

   \item{x.test}{
     Explanatory variables for test (out of sample) data. Should have same column structure as
     \code{x.train}. \code{bart} will generate draws of \eqn{f(x)} for each \eqn{x} which is a
     row of \code{x.test}.
   }

   \item{sigest}{
   	 For continuous response models, an estimate of the error variance, \eqn{\sigma^2},
   	 used to calibrate an inverse-chi-squared prior used on that parameter. If not supplied,
   	 the least-squares estimate is derived instead. See \code{sigquant} for more information. Not applicable when \eqn{y} is binary.
   }

   \item{sigdf}{
     Degrees of freedom for error variance prior. Not applicable when \eqn{y} is binary.
   }

   \item{sigquant}{
     The quantile of the error variance prior that the rough estimate
     (\code{sigest}) is placed at. The closer the quantile is to 1,
     the more aggresive the fit will be as you are putting more prior weight
     on error standard deviations (\eqn{\sigma}) less than the rough estimate.
     Not applicable when \eqn{y} is binary.
   }

   \item{k}{
     For numeric \eqn{y},
     \code{k} is the number of prior standard deviations \eqn{E(Y|x) = f(x)} is away from \eqn{\pm 0.5}{+/- 0.5}.
     The response (\code{y.train}) is internally scaled to range from \eqn{-0.5} to \eqn{0.5}.
     For binary \eqn{y}, \code{k} is the number of prior standard deviations \eqn{f(x)} is away from
     \eqn{\pm 3}{+/- 3}. In both cases, the bigger \eqn{k} is, the more conservative the fitting will be.
   }

   \item{power}{
     Power parameter for tree prior.
   }

   \item{base}{
     Base parameter for tree prior.
   }

   \item{binaryOffset}{
     Used for binary \eqn{y}. When present, the model is
       \eqn{P(Y = 1 \mid x) = \Phi(f(x) + \mathrm{binaryOffset})}{P(Y = 1 | x) = \Phi(f(x) + binaryOffset)},
     allowing fits with probabilities shrunk towards values other than \eqn{0.5}.
   }
   \item{weights}{
      An optional vector of weights to be used in the fitting process. When present, BART
      fits a model with observations \eqn{y \mid x \sim N(f(x), \sigma^2 / w)}{y | x ~ N(f(x), \sigma^2 / w)}, 
      where \eqn{f(x)} is the unknown function.
  }
   \item{ntree}{
     The number of trees in the sum-of-trees formulation.
   }

   \item{ndpost}{
     The number of posterior draws after burn in, \code{ndpost / keepevery} will actually be returned.
   }

   \item{nskip}{
     Number of MCMC iterations to be treated as burn in.
   }

   \item{printevery}{
     As the MCMC runs, a message is printed every \code{printevery} draws.
   }

   \item{keepevery}{
     Every \code{keepevery} draw is kept to be returned to the user. Useful for \dQuote{thinning}
     samples.
   }

   \item{keeptrainfits}{
     If \code{TRUE} the draws of \eqn{f(x)} for \eqn{x} corresponding to the rows of \code{x.train} are returned.
   }

   \item{usequants}{
     When \code{TRUE}, determine tree decision rules using estimated quantiles derived
     from the \code{x.train} variables. When \code{FALSE}, splits are determined using
     values equally spaced across the range of a variable. See details for more information.
   }
   
   \item{numcut}{
   	 The maximum number of possible values used in decision rules (see \code{usequants}, details). If a single
   	 number, it is recycled for all variables; otherwise must be a vector of length equal to
   	 \code{ncol(x.train)}. Fewer rules may be used if a covariate lacks enough unique values.
   }

   \item{printcutoffs}{
     The number of cutoff rules to printed to screen before the MCMC is run.
     Given a single integer, the same value will be used for all variables.
     If 0, nothing is printed.
   }

   \item{verbose}{
     Logical; if \code{FALSE} supress printing.
   }
   
   \item{nchain}{
   	  Integer specifying how many independent tree sets and fits should be calculated.
   }
   
   \item{nthread}{
   	  Integer specifying how many threads to use. Depending on the CPU architecture, using more than the
          number of chains can degrade performance for small/medium data sets. As such
   	  some calculations may be executed single threaded regardless.
   }
   
   \item{keepcall}{
     Logical; if \code{FALSE}, returned object will have \code{call} set to \code{call("NULL")}, otherwise the call
          used to instantiate BART.
   }
   
   \item{formula}{
     The same as \code{x.train}, the name reflecting that a formula object can be used instead.
   }
   
   \item{data}{
     The same as \code{y.train}, the name reflecting that a data frame can be specified when a formula is given instead.
   }
   
   \item{test}{
     The same as \code{x.train}. Can be missing.
   }
   \item{subset}{
     A vector of logicals or indicies used to subset of the data. Can be missing.
   }
   \item{offset}{
     The same as \code{binaryOffset}. Can be missing.
   }
   \item{offset.test}{
     A vector of offsets to be used with test data, in case it is different than the training offset.
     If \code{offest} is missing, defaults to \code{NULL}.
   }
   
   \item{n.trees}{Same as \code{ntree}.}
   \item{n.samples}{Same as \code{ndpost}.}
   \item{n.burn}{Same as \code{nskip}.}
   \item{n.chains}{Same as \code{nchain}.}
   \item{n.threads}{Same as \code{nthread}.}
   \item{n.cuts}{Same as \code{numcut}.}
   \item{n.thin}{Same as \code{keepevery}.}
   \item{keepTrainingFits}{Same as \code{keeptrainfits}.}
   \item{useQuantiles}{Same as \code{usequants}.}
   \item{printEvery}{Same as \code{printevery}.}
   \item{printCutoffs}{Same as \code{printcutoffs}.}
   \item{keepCall}{Same as \code{keepcall}.}

   \item{x}{
     Object of class \code{bart}, returned by function \code{bart}, which contains the information to be plotted.
   }

   \item{plquants}{
     In the plots, beliefs about \eqn{f(x)} are indicated by plotting the
     posterior median and a lower and upper quantile. \code{plquants} is a double vector
     of length two giving the lower and upper quantiles.
   }

   \item{cols}{
     Vector of two colors. First color is used to plot the median of \eqn{f(x)} and
     the second color is used to plot the lower and upper quantiles.
   }

   \item{\dots}{
     Additional arguments passed on to plot or \code{dbartsControl}, respectively.
   }
}
\details{
  BART is an Bayesian MCMC method. 
  At each MCMC interation, we produce a draw from the joint posterior
  \eqn{(f, \sigma) \mid (x, y)}{(f, \sigma) | (x, y)} in the numeric \eqn{y} case
  and just \eqn{f} in the binary \eqn{y} case.

  Thus, unlike a lot of other modeling methods in R, \code{bart} does not produce a single model object
  from which fits and summaries may be extracted. The output consists of values
  \eqn{f^*(x)}{f*(x)} (and \eqn{\sigma^*}{\sigma*} in the numeric case) where * denotes a particular draw.
  The \eqn{x} is either a row from the training data (\code{x.train}) or the test data (\code{x.test}).
  
  \subsection{Decision Rules}{
    Decision rules for any tree are of the form \eqn{x \le c}  vs. \eqn{x > c}
    for each \sQuote{\eqn{x}} corresponding to a column of \code{x.train}.
    \code{usequants} determines the means by which the set of possible \eqn{c} is
    determined. If \code{usequants} is \code{TRUE}, then the \eqn{c} are
    a subset of the values interpolated half-way between the unique, sorted values
    obtained from the corresponding column of \code{x.train}.
    If \code{usequants} is \code{FALSE}, the cutoffs are equally spaced across the
    range of values taken on by the corresponding
    column of \code{x.train}.
     
  	The number of possible values of \eqn{c} is determined by \code{numcut}.
    If \code{usequants} is \code{FALSE}, \code{numcut} equally spaced cutoffs
    are used covering the range of values in the corresponding
    column of \code{x.train}. If \code{usequants} is \code{TRUE}, then for a variable
    the minimum of \code{numcut} and one less than the number
    of unique elements for that variable are used.
  }
}
\value{
  \code{bart} returns a list assigned class \code{bart}. For applicable
  quantities, \code{ndpost / keepevery} samples are returned.
  In the numeric \eqn{y} case, the list has components:

  \item{\code{yhat.train}}{
        A array/matrix of posterior samples. The \eqn{(i, j, k)} value is \eqn{f^*(x_k)}{f*(x_k)}
        for the \eqn{j}th kept draw of of the posterior of \eqn{f}
        and the \eqn{j}th row of \code{x.train} corresponding to chain \eqn{i}. When \code{nchain}
        is one, the result is a collapsed down to a matrix.}
  \item{\code{yhat.test}}{
      Same as \code{yhat.train} but now the \eqn{x}s are the rows of the test data.}
  \item{\code{yhat.train.mean}}{
        Vector of means of \code{yhat.train} across columns and chains, with length equal to the number of
        training observations.}
  \item{\code{yhat.test.mean}}{
        Vector of means of \code{yhat.test} across columns and chains.}
  \item{\code{sigma}}{
        Matrix of posterior samples of \code{sigma}, the residual/error standard deviation. Dimensions
        are equal to the number of chains times the numbers of samples unless \code{nchain} is one.}
  \item{\code{first.sigma}}{
        Burn-in draws of \code{sigma}.}
  \item{\code{varcount}}{
        A matrix with number of rows equal to the number of kept draws and each column
        corresponding to a training variable. Contains the total count of the number of
        times that variable is used in a tree decision rule (over all trees).}
  \item{\code{sigest}}{
        The rough error standard deviation (\eqn{\sigma}) used in the prior.}
  \item{\code{y}}{
        The input dependent vector of values for the dependent variable.
        This is used in \code{plot.bart}.}
    
  In the binary \eqn{y} case, the returned list has the components
  \code{yhat.train}, \code{yhat.test}, and \code{varcount} as above.  In addition the list 
  has a \code{binaryOffset} component giving the value used.

  Note that in the binary \eqn{y}, case \code{yhat.train} and \code{yhat.test} are
  \eqn{f(x) + \mathrm{binaryOffset}}{f(x) + binaryOffset}. For draws of the probability
  \eqn{P(Y = 1 | x)}, apply the normal cdf (\code{pnorm}) to these values.\cr\cr

  The \code{plot} method sets \code{mfrow} to \code{c(1,2)} and makes two plots.
  The first plot is the sequence of kept draws of \eqn{\sigma}
  including the burn-in draws. Initially these draws will decline as BART finds fit
  and then level off when the MCMC has burnt in. The second plot has \eqn{y} on the
  horizontal axis and posterior intervals for the corresponding \eqn{f(x)} on the vertical axis.
}
\references{
Chipman, H., George, E., and McCulloch, R. (2009)
   BART: Bayesian Additive Regression Trees.

Chipman, H., George, E., and McCulloch R. (2006)
   Bayesian Ensemble Learning. 
   Advances in Neural Information Processing Systems 19,
   Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.

both of the above at:
\url{http://www.rob-mcculloch.org}

Friedman, J.H. (1991)
   Multivariate adaptive regression splines.
        \emph{The Annals of Statistics}, \bold{19}, 1--67.
}
\author{
Hugh Chipman: \email{hugh.chipman@gmail.com},\cr
Robert McCulloch: \email{robert.mcculloch1@gmail.com},\cr
Vincent Dorie: \email{vdorie@gmail.com}.
}
\seealso{
\code{\link{pdbart}}
}
\examples{
## simulate data (example from Friedman MARS paper)
## y = f(x) + epsilon , epsilon ~ N(0, sigma)
## x consists of 10 variables, only first 5 matter

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

## run BART
set.seed(99)
bartFit <- bart(x, y)

plot(bartFit)

## compare BART fit to linear matter and truth = Ey
lmFit <- lm(y ~ ., data.frame(x, y))

fitmat <- cbind(y, Ey, lmFit$fitted, bartFit$yhat.train.mean)
colnames(fitmat) <- c('y', 'Ey', 'lm', 'bart')
print(cor(fitmat))
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{nonlinear}
