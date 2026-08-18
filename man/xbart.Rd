\name{xbart}
\alias{xbart}
\title{Crossvalidation For Bayesian Additive Regression Trees}
\description{
  Fits the BART model against varying \code{k}, \code{power}, \code{base}, and \code{n.trees} parameters using \eqn{K}-fold or repeated random subsampling crossvalidation, sharing burn-in between parameter settings. Results are returned as an array of evaluations of a loss function on the held-out sets.
}
\usage{
xbart(
    formula, data, subset, weights, offset, verbose = FALSE, n.samples = 200L,
    method = c("k-fold", "random subsample"), n.test = c(5, 0.2),
    n.reps = 40L, n.burn = c(200L, 150L),
    loss = c("rmse", "log", "mcr"), n.threads = dbarts::guessNumCores(), n.trees = 75L,
    k = NULL, power = 2, base = 0.95,
    split.probs = NULL, dart = FALSE, drop = TRUE,
    resid.prior = chisq, sigest = NA_real_,
    seed = NA_integer_,
    factors = c("categorical", "indicators"),
    family = c("auto", "gaussian", "probit", "logistic"),
    missing = c("incorporate", "error"),
    node.prior = NULL, n.cuts = 100L, useQuantiles = FALSE, n.thin = 1L,
    storage = c("double", "single"), tree.prior = NULL)
}
\arguments{
  \item{formula}{
    An object of class \code{\link{formula}} following an analogous model description syntax as \code{\link{lm}}. For backwards compatibility, can also be the \code{\link{bart}} matrix \code{x.train}. See \code{\link{dbarts}}.
  }
  \item{data}{
    An optional data frame, list, or environment containing predictors to be used with the model. For backwards compatibility, can also be the \code{\link{bart}} vector \code{y.train}.
  }
  \item{subset}{
    An optional vector specifying a subset of observations to be used in the fitting process.
  }
  \item{weights}{
    An optional vector of weights to be used in the fitting process. For a gaussian response, BART fits a model with observations \eqn{y \mid x \sim N(f(x), \sigma^2 / w)}{y | x ~ N(f(x), \sigma^2 / w)}, where \eqn{f(x)} is the unknown function. Binary responses differ: a \code{"probit"} model does not support weights (a weighted probit has no tractable latent-variable form), except that weights identically 1 are treated as absent; a \code{"logistic"} model treats them as observation counts and so requires positive integers (its Polya-Gamma latent for a count \eqn{w} is a sum of \eqn{w} unit draws).
  }
  \item{offset}{
    An optional vector specifying an offset from 0 for the relationship between the underlying function, \eqn{f(x)}, and the response \eqn{y}. Useful for all three response families, though each applies it differently: for a gaussian response, \eqn{y = f(x) + \mathrm{offset} + \epsilon}{y = f(x) + offset + epsilon}, a fixed component of the mean (see \code{\link{dbarts}} for the interaction with BART's internal range-scaling); for binary responses it enters the link directly, \eqn{P(Y = 1 \mid X = x) = \Phi(f(x) + \mathrm{offset})}{P(Y = 1 | X = x) = \Phi(f(x) + offset)} for a \code{"probit"} fit and the analogous model on the logistic link for \code{"logistic"}.
  }
  \item{verbose}{
    A logical determining if additional output is printed to the console.
  }
  \item{n.samples}{
    A positive integer, setting the number of posterior samples drawn for each fit of training data and used by the loss function.
  }
  \item{method}{
     Character string, either \code{"k-fold"} or \code{"random subsample"}.
   }
  \item{n.test}{
    For each fit, the test sample size or proportion. For method \code{"k-fold"}, is expected to be the number of folds, and in \eqn{[2, n]}. For method \code{"random subsample"}, can be a real number in \eqn{(0, 1)} or a positive integer in \eqn{(1, n)}. When a given as proportion, the number of test observations used is the proportion times the sample size rounded to the nearest integer.
  }
  \item{n.reps}{
    A positive integer setting the number of cross validation steps that will be taken. For \code{"k-fold"}, each replication corresponds to fitting each of the \eqn{K} folds in turn, while for \code{"random subsample"} a replication is a single fit.
  }
  \item{n.burn}{
    Unlike the single-scalar \code{n.burn} of \code{\link{dbarts}}, \code{bart2}, and \code{rbart_vi}, here it is one or two non-negative integers, specifying 1) the burn-in when a chain is freshly started against a data split and 2) the burn-in when moving from one parameter setting to another over the same split. Chains are never carried between data splits or folds - the held-out observations of one were training observations of the previous, so continuing a chain lets slowly-mixing settings score against data they have effectively seen.
  }
  \item{loss}{
    Either one of the pre-set loss functions as character-strings (\code{mcr} - misclassification rate for binary responses, \code{rmse} - root-mean-squared-error for continuous response), \code{log} - negative log-loss for binary response (\code{rmse} serves this purpose for continuous responses), a function, or a function-evaluation environment list-pair. Functions should have prototypes of the form \code{function(y.test, y.test.hat, weights)}, where \code{y.test} is the held out test subsample, \code{y.test.hat} is a matrix of dimension \code{length(y.test)} \eqn{\times}{*} \code{n.samples}, and \code{weights} are an optional vector of user-supplied weights. See examples.
  }
  \item{n.threads}{
    Replications are independent, and for \code{n.threads > 1} they are divided into approximately equal chunks and executed on that many parallel workers (a \code{\link[parallel]{makeCluster}} cluster). The default uses \code{\link{guessNumCores}}, which should work across the most common operating system/hardware pairs.
   }
   \item{n.trees}{
     A vector of positive integers setting the BART hyperparameter for the number of trees in the sum-of-trees formulation. See \code{\link{bart}}.
   }
   \item{k}{
     A vector of positive real numbers, setting the BART hyperparameter for the node-mean prior standard deviation. If \code{NULL}, the default of \code{bart2} will be used - 2 for continuous response and a Chi hyperprior for binary. Hyperprior crossvalidation not possible at this time. A hyperprior \code{k} is held, not swept, and is DRAWN every sweep in every cell, so the reported loss is computed under a shrinkage that moves within each fit rather than under the named value.
   }
   \item{power}{
     A vector of real numbers greater than one, setting the BART hyperparameter for the tree prior's growth probability, given by \eqn{{base} / (1 + depth)^{power}}.
   }
   \item{base}{
     A vector of real numbers in \eqn{(0, 1)}, setting the BART hyperparameter for the tree prior's growth probability.
   }
   \item{split.probs}{
     Prior probabilities that a variable is used in a splitting rule, as in \code{\link{bart2}}. A single value or \code{NULL} yields the uniform default; a named or unnamed vector assigns per-column probabilities. Fixed for the whole crossvalidation, not part of the swept grid. Cannot be combined with \code{dart}.
   }
   \item{dart}{
     Use the DART sparsity-inducing prior of Linero (2018) instead of a fixed \code{split.probs}. \code{TRUE} enables it with default hyperparameters; a prior created by \code{dbartsPriors$dart} supplies its own. The split probabilities are then sampled rather than fixed, while \code{power} and \code{base} are still swept.
   }
   \item{drop}{
     Logical, determining if dimensions with a single value are dropped from the result.
   }
   \item{resid.prior}{
     An expression of the form \code{chisq} or \code{chisq(df, quant)} that sets the prior used on the residual/error variance.
   }
   \item{sigest}{
     A positive numeric estimate of the residual standard deviation. If \code{NA}, a linear model is used with all of the predictors to obtain one. Fitting functions (\code{xbart}, \code{\link{bart}}/\code{bart2}, \code{\link{rbart_vi}}) spell this \code{sigest}; sampler constructors (\code{\link{dbarts}}, \code{dbartsSpec}) spell the same concept \code{sigma}.
   }
   \item{seed}{
     Optional integer specifying the desired pRNG \link[=set.seed]{seed}. Results are reproducible for a fixed pair of \code{seed} and \code{n.threads}; the caller's random stream is left untouched when a seed is given. Without one, \code{\link{set.seed}} beforehand suffices.
   }
   \item{factors}{
     How factor columns in a data frame enter the model; as in \code{\link{dbarts}}.
   }
   \item{family}{
     The response model, resolved as in \code{\link{dbarts}}: \code{"auto"} fits gaussian models to continuous responses and probit models to those coded 0/1, \code{"gaussian"} forces a continuous fit, and \code{"probit"} and \code{"logistic"} require a 0/1 response. A two-level factor, logical, or two-level character response is detected and fit as probit; a factor with three or more levels is an error, as \code{xbart} does not cross-validate the multinomial model. The built-in binary losses transform test predictions through the family's link.
   }
   \item{missing}{
     How missing values in the predictors enter the model; as in \code{\link{dbarts}}.
   }
   \item{node.prior}{
     An optional expression of the form \code{normal(k)}, \code{linear(columns, k)}, or \code{gp(columns, k, ...)} selecting the leaf model, as in \code{\link{dbarts}}. The default fits constant leaves. A \code{k} given inside the prior stands in for a missing \code{k} argument; the \code{k} argument otherwise drives the crossvalidation grid as usual.

     A named leaf calibration (\code{normal(scale = )}, and likewise for \code{linear} and \code{gp}) is NOT a grid axis: it is held fixed across every cell, so the \code{k} grid sweeps the prior standard deviation \code{scale / k} about a fixed anchor rather than about whatever the response range happened to imply. It is honored in cells that re-use a sampler as well as in cells that create one, so the loss surface does not depend on the order the cells run in. The \code{sd} spelling meets the same refusal here as everywhere when \code{k} is a hyperprior.
   }
   \item{n.cuts}{
     A positive integer giving the number of decision rules used for each predictor, as in \code{\link{dbartsControl}}.
   }
   \item{useQuantiles}{
     Logical; as in \code{\link{dbartsControl}}, determines whether decision rules are placed at the empirical quantiles of each predictor's values rather than spaced uniformly through its range.
   }
   \item{n.thin}{
     A positive integer; as in \code{\link{dbartsControl}}, thins each cell's chain against serial correlation. \code{n.samples} are still returned regardless.
   }
   \item{storage}{
     A character string selecting the precision of the internal running residual, spelled and defaulted as in \code{\link{dbartsControl}}. Every cell's sampler is created over a shared per-fold data handle, a path the engine keeps in double precision regardless of family or leaf model; \code{"single"} is refused rather than silently ignored.
   }
   \item{tree.prior}{
     An expression of the form \code{cgm} or \code{cgm(power, base, split.probs)} that sets the tree structure prior, or a prior object built with \code{\link{dbartsPriors}} - \code{cgm(...)} or \code{dbartsPriors$dart(...)}. \code{power} and \code{base} are xbart's grid axes, so a supplied object's own \code{power}/\code{base} are replaced by the swept grid values every cell exactly as the \code{k} argument replaces a supplied \code{node.prior}'s \code{k}; the object's other content - a \code{cgm} object's \code{split.probs}, a DART object's Dirichlet hyperparameters - rides every cell unchanged. Because \code{power}/\code{base} are grid axes here rather than ordinary scalars, they may be supplied alongside \code{tree.prior}; \code{dart} and \code{split.probs} would only duplicate what a supplied \code{tree.prior} already specifies, so combining either with it is an error naming both.
   }
}
\details{
  Crossvalidates \code{n.reps} replications against the crossproduct of given hyperparameter vectors \code{n.trees} \eqn{\times}{*} \code{k} \eqn{\times}{*} \code{power} \eqn{\times}{*} \code{base}. For each fit, either one fold is withheld as test data and \code{n.test - 1} folds are used as training data or \code{n * n.test} observations are withheld as test data and \code{n * (1 - n.test)} used as training. A replication corresponds to fitting all \eqn{K} folds in \code{"k-fold"} crossvalidation or a single fit with \code{"random subsample"}. The training data is used to fit a model and make predictions on the test data which are used together with the test data itself to evaluate the \code{loss} function.

  \code{loss} functions are either the default of average negative log-loss for binary outcomes and root-mean-squared error for continuous outcomes, misclassification rates for binary outcomes, or a \code{function} with arguments \code{y.test} and \code{y.test.hat}. \code{y.test.hat} is of dimensions equal to \code{length(y.test)} \eqn{\times}{*} \code{n.samples}, on the latent scale for binary outcomes. A third option is to pass a list of \code{list(function, evaluationEnvironment)}, so as to provide default bindings. The binary losses \code{"log"} and \code{"mcr"} cannot apply to continuous responses.

}
\value{
  An array of dimensions \code{n.reps} \eqn{\times}{*} \code{length(n.trees)} \eqn{\times}{*} \code{length(k)} \eqn{\times}{*} \code{length(power)} \eqn{\times}{*} \code{length(base)}. If \code{drop} is \code{TRUE}, dimensions of length 1 are omitted. If all hyperparameters are of length 1, then the result will be a vector of length \code{n.reps}. When the result is an array, the \code{dimnames} of the result shall be set to the corresponding hyperparameters.

  For method \code{"k-fold"}, each element is an average across the \eqn{K} fits. For \code{"random subsample"}, each element represents a single fit.
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

mad <- function(y.test, y.test.hat, weights) {
    # note, weights are ignored
    mean(abs(y.test - apply(y.test.hat, 1L, mean)))
}



## low iteration numbers to to run quickly
xval <- xbart(x, y, n.samples = 15L, n.reps = 4L, n.burn = c(10L, 3L),
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

