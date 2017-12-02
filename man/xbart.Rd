\name{xbart}
\alias{xbart}
\title{Crossvalidation For Bayesian Additive Regression Trees}
\description{
Fits the BART model against varying \code{k}, \code{power}, \code{base}, and \code{ntree} parameters using \eqn{K}-fold crossvalidation, sharing burn-in between parameter settings. Results are given an array of evalulations of a loss functions on the held-out sets.
}
\usage{
xbart(formula, data, subset, weights, offset, verbose = FALSE,
      n.samples = 200L, K = 5L, n.reps = 200L, n.burn = c(200L, 150L, 50L),
      loss = c("rmse", "mcr"), n.threads = guessNumCores(), n.trees = 75L,
      k = 2, power = 2, base = 0.95, drop = TRUE,
      resid.prior = chisq, control = dbartsControl(), sigma = NA_real_)}
\arguments{
  \item{formula}{
    An object of class \code{\link{formula}} following an analogous model description
    syntax as \code{\link{lm}}. For backwards compatibility, can also be the \code{\link{bart}} matrix
    \code{x.train}. See \code{\link{dbarts}}.
  }
  \item{data}{
    An optional data frame, list, or environment containing predictors to be used with the
    model. For backwards compatibility, can also be the \code{\link{bart}} vector \code{y.train}.
  }
  \item{subset}{
    An optional vector specifying a subset of observations to be used in the fitting process.
  }
  \item{weights}{
    An optional vector of weights to be used in the fitting process. When present, BART
    fits a model with observations \eqn{y \mid x \sim N(f(x), \sigma^2 / w)}{y | x ~ N(f(x), \sigma^2 / w)}, 
    where \eqn{f(x)} is the unknown function.
  }
  \item{offset}{
    An optional vector specifying an offset from 0 for the relationship between the underyling function,
    \eqn{f(x)}, and the response \eqn{y}. Only is useful for binary responses, in which case the model
    fit is to assume
    \eqn{P(Y = 1 \mid X = x) = \Phi(f(x) + \mathrm{offset})}{P(Y = 1 | X = x) = \Phi(f(x) + offset)},
    where \eqn{\Phi} is the standard normal cumulative distribution function.
  }
  \item{verbose}{
    A logical determining if additional output is printed to the console.
  }
  \item{n.samples}{
    A positive integer, setting the number of posterior samples drawn for each replication and used by
    the loss function.
  }
  \item{K}{
    Number of \eqn{K}-folds for use in crossvalidation. In each replication, \eqn{1/K} of the data will be
    held out and used for evaluation, while \eqn{1 - 1/K} are used for training.
  }
  \item{n.reps}{
    A positive integer setting the number of \eqn{K}-fold crossvalidation evaluations of the loss function
    will be evaluated.
  }
  \item{n.burn}{
    Between one and three positive integers, specifying the 1) initial burn-in, 2) burn-in when moving from
    one parameter setting to another, and 3) the burn-in between each \eqn{K}-fold replication.
  }
  \item{loss}{
    Either a one of the present loss functions as character-strings (\code{mcr} - missclassification rate for
    binary responses, \code{rmse} - root-mean-squared-error for continuous response), a function, or a function-
    evaluation environment list-pair. Functions should have prototypes of the form
    \code{function(y.test, y.test.hat)}, where \code{y.test} is the held out \eqn{K}-fold and \code{y.test.hat}
    is a matrix of dimension \code{length(y.test) x n.samples}. See examples.
  }
  \item{n.threads}{
    Across different sets of parameters (\code{k} \eqn{\times}{*} \code{power} \eqn{\times}{*} \code{base}
    \eqn{\times}{*} \code{n.trees}) and \code{n.reps}, results are independent. For \code{n.threads > 1},
    evaluations of the above are divided into approximately equal size evaluations chunks and executed in parallel.
    The default uses \code{link{guessNumCores}}, which should work across the most common operating
    system/hardware pairs. A value of \code{NA} is interpetted as 1.
   }
   \item{n.trees}{
     A vector of positive integers setting the BART hyperparameter for the number of trees in the
     sum-of-trees formulation. See \code{\link{bart}}.
   }
   \item{k}{
     A vector of positive real numbers, setting the BART hyperparameter for the node-mean prior standard
     deviation.
   }
   \item{power}{
     A vector of real numbers greater than one, setting the BART hyperparameter for the tree prior's growth
     probability, given by \eqn{{base} / (1 + depth)^{power}}.
   }
   \item{base}{
     A vector of real numbers in \eqn{(0, 1)}, setting the BART hyperparameter for the tree prior's growth
     probability.
   }
   \item{drop}{
     Logical, determining if dimensions with a single value are dropped from the result.
   }
   \item{resid.prior}{
     An expression of the form \code{chisq} or \code{chisq(df, quant)} that sets the prior used on
     the residual/error variance.
   }
   \item{control}{
     An object inheriting from \code{dbartsControl}, created by the \code{\link{dbartsControl}} function.
   }
   \item{sigma}{
     A positive numeric estimate of the residual standard deviation. If \code{NA}, a linear model is used
     with all of the predictors to obtain one.
   }
}
\details{
  Fits \eqn{K}-fold cross validation for \code{n.reps} replications against the crossproduct of given
  hyperparameter vectors \code{n.trees} \eqn{\times}{*} \code{k} \eqn{\times}{*} \code{power}
  \eqn{\times}{*} \code{base}. \eqn{1/K} of the observations are withheld from each replication while
  \eqn{1 - 1/K} are used to train the BART model, and for each the withheld observations are used
  together with BART estimates to calculate the \code{loss} function.
  
  \code{loss} functions are either the default of missclassification rates for binary outcomes
  or a \code{function} with arguments \code{y.test}, \code{y.test.hat} that returns a scalar value.
  \code{y.test.hat} is of dimensions equal to \code{length(y.test)} \eqn{\times}{*} \code{n.samples}.
  A third option is to pass a list of \code{list(function, evaluationEnvironment)}, so as to
  provide default bindings.
}
\value{
  An array of dimensions \code{n.reps} \eqn{\times}{*} \code{length(n.trees)} \eqn{\times}{*}
  \code{length(k)} \eqn{\times}{*} \code{length(power)} \eqn{\times}{*} \code{length(base)}.
  If \code{drop} is \code{TRUE}, dimensions of length 1 are omitted. If all hyperparameters
  are of length 1, then the result will be a vector of length \code{n.reps}. When the result is an
  array, the \code{dimnames} of the result shall be set to the corresponding hyperparameters.
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

mad <- function(y.train, y.train.hat) 
    mean(abs(y.train - apply(y.train.hat, 1L, mean)))



## low iteration numbers to to run quickly
xval <- xbart(x, y, n.samples = 15L, n.reps = 8L, n.burn = c(10L, 3L, 1L),
              n.trees = c(5L, 7L),
              k = c(1, 2, 4),
              power = c(1.5, 2),
              base = c(0.75, 0.8, 0.95), n.threads = 1L,
              loss = mad)
}
\keyword{nonparametric}
\keyword{tree}
\keyword{regression}
\keyword{crossvalidation}
