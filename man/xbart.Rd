\name{xbart}
\alias{xbart}
\title{Crossvalidation For Bayesian Additive Regression Trees}
\description{
  Fits the BART model against varying \code{k} (or \code{sd}), \code{power}, \code{base}, and \code{n.trees} parameters using \eqn{K}-fold or repeated random subsampling crossvalidation, sharing burn-in between parameter settings. Results are returned as an array of evaluations of a loss function on the held-out sets.
}
\usage{
xbart(
    formula, data, subset, weights, offset, verbose = FALSE, n.samples = 200L,
    method = c("k-fold", "random subsample"), n.test = c(5, 0.2),
    n.reps = 40L, n.burn = c(200L, 150L),
    loss = c("rmse", "log", "mcr"), n.threads = dbarts::guessNumCores(), n.trees = 75L,
    k = NULL, sd = NULL, power = 2, base = 0.95,
    split.probs = NULL, drop = TRUE,
    sigest = NULL,
    seed = NULL,
    factors = c("categorical", "indicators"),
    family = c("auto", "gaussian", "probit", "logistic"),
    leaf.prior = NULL, n.cuts = 100L, useQuantiles = FALSE, n.thin = 1L,
    storage = c("double", "single"), tree.prior = NULL,
    parallel = getOption("dbarts.parallel", "auto"), cl = NULL,
    control = dbarts::dbartsControl(), sigma = NULL, \dots)
}
\arguments{
  \item{sigma}{
    The 0.9-x spelling of \code{sigest}, accepted for one release with a once-per-session warning and removed in dbarts 1.1-0. Supplying both is an error.
  }
  \item{\dots}{
    Not used for new code: the channel that lets a retired argument spelling (\code{dart}, now \code{tree.prior = dart()}; \code{resid.prior}, now \code{family = gaussian(sigma = )}, its one home, so writing it both ways is refused where the two disagree) reach a message naming its successor, instead of R's own \dQuote{unused argument} error. Any other name is refused. Removed in dbarts 1.1-0.
  }
  \item{control}{
    A \code{\link{dbartsControl}} object carrying the sampler and engine settings, including the ones \code{xbart} spells no flat name for: \code{categoricalExhaustiveCap}, \code{testFitParallelCutoff}, \code{predictParallelCutoff}, \code{sparseDensityThreshold}, \code{treeShift}, and the tree-move mixture \code{proposal.probs}. Precedence, one rule: a flat argument named in the call wins over the control's slot of the same name, and a slot the control speaks for - one its own \code{dbartsControl()} call named, or one edited afterwards to differ from a fresh control's - wins over this function's default, while a slot it never spoke for leaves that default standing. The fields a sweep forces - \code{n.chains}, the control's own \code{n.threads}, \code{keepTrees}, \code{keepTrainingFits}, \code{updateState}, \code{verbose} - are applied after both, and \code{n.burn} and \code{n.threads} here are the sweep's own grid axis and worker count rather than the control fields of the same name. A control taken from a fitted sampler carries that fit's model configuration and is refused by name - pass a fresh \code{dbartsControl()}.
  }
  \item{formula}{
    An object of class \code{\link{formula}} following an analogous model description syntax as \code{\link{lm}}. For backwards compatibility, can also be the \code{\link{bartBT}} matrix \code{x.train}. See \code{\link{dbarts}}. A \code{dbartsData} object carrying \code{bases} (a multi-forest model) is refused: \code{xbart} cross-validates single-forest models.
  }
  \item{data}{
    An optional data frame, list, or environment containing predictors to be used with the model. For backwards compatibility, can also be the \code{\link{bartBT}} vector \code{y.train}.
  }
  \item{subset}{
    An optional vector specifying a subset of observations to be used in the fitting process.
  }
  \item{weights}{
    An optional vector of weights to be used in the fitting process. For a gaussian response, BART fits a model with observations \eqn{y \mid x \sim N(f(x), \sigma^2 / w)}{y | x ~ N(f(x), \sigma^2 / w)}, where \eqn{f(x)} is the unknown function. Binary responses differ: a \code{"probit"} model does not support weights (a weighted probit has no tractable latent-variable form), except that weights identically 1 are treated as absent - the all-0/1 vector \code{\link{dbarts}} installs as an active-row mask is refused here, cross-validation partitioning the rows itself; a \code{"logistic"} model treats them as observation counts and so requires positive integers (its Polya-Gamma latent for a count \eqn{w} is a sum of \eqn{w} unit draws, so a sweep's time grows with the total count, and a count above \eqn{10^6}{1e6} is refused).
  }
  \item{offset}{
    An optional vector specifying an offset from 0 for the relationship between the underlying function, \eqn{f(x)}, and the response \eqn{y}. Useful for all three response families, though each applies it differently: for a gaussian response, \eqn{y = f(x) + \mathrm{offset} + \epsilon}{y = f(x) + offset + epsilon}, a fixed component of the mean (see \code{\link{dbarts}} for the interaction with BART's internal range-scaling); for binary responses it enters the link directly, \eqn{P(Y = 1 \mid X = x) = \Phi(f(x) + \mathrm{offset})}{P(Y = 1 | X = x) = \Phi(f(x) + offset)} for a \code{"probit"} fit and the analogous model on the logistic link for \code{"logistic"}.
  }
  \item{verbose}{
    A logical determining if additional output is printed to the console, including the one-line message \code{family = "auto"} prints, once per call, naming the family it resolves to.
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
    Unlike the single-scalar \code{n.burn} of \code{\link{dbarts}} and \code{bart}, here it is one or two non-negative integers, specifying 1) the burn-in when a chain is freshly started against a data split and 2) the burn-in when moving from one parameter setting to another over the same split. A longer vector is an error naming the argument rather than being silently truncated to its first two entries; dbarts 0.9-x read a third element as a per-replication burn-in, which no longer exists. Chains are never carried between data splits or folds - the held-out observations of one were training observations of the previous, so continuing a chain lets slowly-mixing settings score against data they have effectively seen.
  }
  \item{loss}{
    Either one of the pre-set loss functions as character-strings (the default is \code{"rmse"} for a continuous response and \code{"log"} for a binary one; \code{mcr} - misclassification rate for binary responses, \code{rmse} - root-mean-squared-error for continuous response), \code{log} - negative log-loss for binary response (\code{rmse} serves this purpose for continuous responses), a function, or a function-evaluation environment list-pair. A function keeps its own environment, so a closure reads what it captured; the list form calls the function from the given environment. Functions should have prototypes of the form \code{function(y.test, y.test.hat, weights)}, where \code{y.test} is the held out test subsample, \code{y.test.hat} is a matrix of dimension \code{length(y.test)} \eqn{\times}{*} \code{n.samples}, and \code{weights} are an optional vector of user-supplied weights. See examples.
  }
  \item{n.threads}{
    Every (replication, fold) pair is an independent unit of work, and for \code{n.threads > 1} the units are divided into approximately equal chunks and executed on that many parallel workers (see \code{parallel}). A \code{k}-fold run of a single replication therefore uses up to \code{k} workers. The default uses \code{\link{guessNumCores}}, which should work across the most common operating system/hardware pairs, and is one where the cores cannot be counted; a stated \code{NA} is refused. Warnings raised while fitting or scoring, including by a supplied \code{loss} function, are collected from every unit and signalled in the calling session once all units have finished, in unit order and each distinct warning once, so the same warnings are seen at any \code{n.threads}.
   }
   \item{parallel}{
    One of \code{"auto"}, \code{"fork"} or \code{"socket"}, and by default \code{getOption("dbarts.parallel", "auto")}. \code{"fork"} runs the workers as forked copies of the calling session, which is fast to start and shares the data; \code{"socket"} starts separate R sessions (a \code{\link[parallel]{makeCluster}} cluster) that each load dbarts and receive their own copy of the data, which costs roughly 80 MB per worker beyond the data itself. \code{"auto"} forks unless the platform is Windows, the session runs under RStudio or Positron, or R.app, in which case it uses sockets. \code{"fork"} is an error on Windows. Forking can hang or crash when the calling session has already used a library that is not fork-safe, such as a multithreaded BLAS; if so, set \code{parallel = "socket"} or \code{options(dbarts.parallel = "socket")}. Results do not depend on the kind of worker. Ignored when \code{cl} is given or only one worker is used.
  }
  \item{cl}{
    \code{NULL} or a cluster from the \code{parallel} package that the caller has already started, used instead of starting workers whatever \code{parallel} says. The units are divided into \code{min(n.threads, length(cl))} chunks, and \code{xbart} does not stop the cluster. The workers must load the same dbarts version as the calling session, since the work function is sent by reference and each worker runs its own installed dbarts. A one-node \code{cl}, like \code{n.threads = 1}, runs the units in the calling session.
  }
  \item{n.trees}{
     A vector of positive integers setting the BART hyperparameter for the number of trees in the sum-of-trees formulation. See \code{\link{bart}}.
   }
   \item{k}{
     The grid for the BART hyperparameter setting the leaf-mean prior standard deviation: a vector of positive real numbers, each a cell of its own, or a \code{\link{list}} whose entries are positive numbers and \code{k} hyperpriors (\code{\link{dbartsPriors}}), so that \code{k = list(1, 2, chi())} sweeps two fixed values against the modelled one. A modelled cell holds its hyperprior and DRAWS \code{k} every sweep, so its loss is computed under a shrinkage that moves within the fit rather than under a named value; it is labelled in the result by the constructor call that rebuilds it. If \code{NULL}, one cell is run at the default \code{\link{bart}} would fit for the response type: fixed 2 for a continuous response and \code{chi(1.5, 2)} for a binary one. A \code{k} carried by a supplied \code{leaf.prior} stands in for a missing argument. Fixed cells are always swept largest to smallest, with a modelled cell last, and the reported \code{k} axis is un-permuted back to the order given, so results do not depend on the order \code{k} is listed in; cells still warm-start off the previous one, so this is order-invariance, not an unbiased estimate for each cell taken alone.
   }
   \item{sd}{
     The same grid axis stated in absolute spreads instead, exclusive with \code{k}: a vector of leaf-prior sds on the family's scale (see \code{normal(sd = )} in \code{\link{dbartsPriors}}), or a \code{\link{list}} mixing them with \code{invchi} laws. Each spread is held fixed across folds, where a \code{k} is relative to each fold's own training range; on a binary family, whose \code{k.scale} is a constant of the latent scale, a grid of fixed cells is the same sweep either way (\code{sd = k.scale / k}); a modelled cell is not, since \code{invchi(df, c)} starts its chain at the spread \code{c} where \code{chi(df, k.scale / c)} starts at \code{k.scale / 2}. A \code{k} inside \code{leaf.prior} beside an \code{sd} grid is refused, as an \code{sd} there is beside either grid. The cells are swept most shrunk first - the smallest sd - with the warm starts and unit parallelism the \code{k} grid has, and the result's axis is labelled \code{sd}.
   }
   \item{power}{
     A vector of real numbers greater than one, setting the BART hyperparameter for the tree prior's growth probability, given by \eqn{{base} / (1 + depth)^{power}}.
   }
   \item{base}{
     A vector of real numbers in \eqn{(0, 1)}, setting the BART hyperparameter for the tree prior's growth probability.
   }
   \item{split.probs}{
     Prior probabilities that a variable is used in a splitting rule, as in \code{\link{bart}}. A single value or \code{NULL} yields the uniform default; a named or unnamed vector assigns per-column probabilities. Fixed for the whole crossvalidation, not part of the swept grid. Cannot be combined with a DART \code{tree.prior}.
   }
   \item{drop}{
     Logical, determining if dimensions with a single value are dropped from the result.
   }
   \item{sigest}{
     A positive numeric estimate of the residual standard deviation. If \code{NULL} (the default), a linear model is used with all of the predictors to obtain one, fit separately for each fold on that fold's training rows, so no fold's prior reads its held-out responses (a fit on all rows runs once beforehand to raise any refusal or fallback once; where that fit falls back to the marginal standard deviation, as described below, each fold takes its own training rows' marginal standard deviation without repeating the linear model); an explicit \code{NA} is a missing value and is refused. Every entry point spells this \code{sigest} - \code{xbart}, \code{\link{bartBT}}/\code{bart} and the sampler constructors \code{\link{dbarts}}/\code{dbartsSpec} alike; \code{sigma} is the retired 0.9-x spelling on \code{dbarts} and on \code{xbart} itself (\code{dbartsSpec} refuses it); a family object's own \code{sigma} is the prior this estimate calibrates, not the estimate. That estimate falls back to the marginal standard deviation of the response when the linear model's residual standard error comes out non-finite, warning as it does so (class \code{dbartsSigmaFallbackWarning}); a design with sparse-backed predictor columns skips the linear model altogether and falls back the same way (class \code{dbartsSparseSigmaFallbackWarning}, a \code{dbartsSigmaFallbackWarning}). It is the estimate a \code{chisq} residual prior's quantile is calibrated against, so it stands beside \code{family = gaussian(sigma = chisq(df, quant))} and is refused beside \code{family = gaussian(sigma = fixed(value))}, which fixes the residual scale outright and would overwrite the estimate with its square root.
   }
   \item{seed}{
     Optional integer specifying the desired pRNG \link[=set.seed]{seed}. \code{NULL} (the default) means not given here and defers to a seed already sitting in \code{control}, if any; an \code{NA} is a missing value, as in \code{\link{set.seed}}: it reads as \code{NULL} for one release, with a once-per-session warning. From the seed in force, \code{xbart} draws a split seed for each replication and a seed for each sampler its units create - one per (replication, fold) unit at each distinct tree count - with \code{\link{sample.int}}, in one pass, under the caller's own \code{\link{RNGkind}}: a given seed therefore draws different values, and so gives different results, under a different kind. Each unit's sampler then takes its seed through \code{control}'s seed slot, exactly as \code{\link{bart}}'s own \code{seed} does, driving a dedicated generator that never reads R's stream, so results are reproducible at any \code{n.threads} for a given seed and \code{RNGkind}. A supplied \code{loss} function that draws random numbers draws on the worker running it, from that process's own default generator; at \code{n.threads} greater than 1 this makes the run not reproducible, whatever seed \code{xbart} itself is given. The caller's random stream is left untouched when a seed is given, and advanced only by that derivation when one is not; without a seed, \code{\link{set.seed}} beforehand suffices. See the Reproducibility section of \link{bart}.
   }
   \item{factors}{
     How factor columns in a data frame enter the model; as in \code{\link{dbarts}}.
   }
   \item{family}{
     The response model, resolved as in \code{\link{dbarts}}: \code{"auto"} fits gaussian models to continuous responses and probit models to those coded 0/1, \code{"gaussian"} forces a continuous fit, and \code{"probit"} and \code{"logistic"} require a 0/1 response. A two-level factor, logical, or two-level character response is detected and fit as probit; a factor with three or more levels is an error, as \code{xbart} does not cross-validate the multinomial model, and neither is a matrix response (per-category counts or a (time, status) pair), which is refused naming \code{bart}. A survival response - a \code{Surv} object or two-column \code{(time, status)} pair, on \code{formula} or as \code{data} directly - is likewise an error naming \code{bart}/\code{dbarts} with \code{family = "aft"}/\code{"hazard"}: \code{xbart} does not cross-validate a survival model. The built-in binary losses transform test predictions through the family's link. This vocabulary is narrower than \code{\link{bart}}'s by design - \code{xbart} cross-validates a single scalar loss per fold, which the own-class families' K-forest or two-part fits have no single counterpart of; the wider family set lives on \code{bart}. Base R family objects map as \code{\link{glm}} takes them (\code{binomial} is the logit link, so \code{family = binomial} is \code{"logistic"}, not probit) and the rest are refused; see \code{\link{dbartsFamilies}}.
   }
   \item{leaf.prior}{
     An optional expression of the form \code{normal(k)}, \code{linear(columns, k)}, or \code{gp(columns, k, ...)} selecting the leaf model, as in \code{\link{dbarts}}. The default fits constant leaves. A \code{k} given inside the prior stands in for a missing \code{k} argument; the \code{k} argument otherwise drives the crossvalidation grid as usual.

     A named \code{sd} inside the prior (\code{normal(sd = )}, and likewise for \code{linear} and \code{gp}) stands in for a missing \code{sd} argument, as a one-cell \code{sd} axis. Beside a \code{k} or \code{sd} grid it is refused, since both would state the spread: name the spreads in the \code{sd} grid. Two calls with the same \code{seed} use the same folds whatever their spreads, so one call per \code{sd} reproduces each cell a grid's fresh start would.
   }
   \item{n.cuts}{
     A positive integer giving the number of decision rules used for each predictor, as in \code{\link{dbartsControl}}.
   }
   \item{useQuantiles}{
     Logical; as in \code{\link{dbartsControl}}, determines whether decision rules are placed at the empirical quantiles of each predictor's values rather than spaced uniformly through its range. Either way a fold takes the full data's rules and their weights, so on a predictor with fewer distinct values than \code{n.cuts} it holds a rule between held-out values as well.
   }
   \item{n.thin}{
     A positive integer; as in \code{\link{dbartsControl}}, thins each cell's chain against serial correlation. \code{n.samples} are still returned regardless.
   }
   \item{storage}{
     A character string selecting the precision of the internal running residual, spelled and defaulted as in \code{\link{dbartsControl}}. Every cell's sampler is created over a shared per-fold data handle, a path the engine keeps in double precision regardless of family or leaf model; \code{"single"} is refused rather than silently ignored.
   }
   \item{tree.prior}{
     An expression of the form \code{cgm} or \code{cgm(power, base, split.probs)} that sets the tree structure prior, or a prior object built with \code{\link{dbartsPriors}} - \code{cgm(...)} or \code{dbartsPriors$dart(...)}. \code{power} and \code{base} are xbart's grid axes, so a supplied object's own \code{power}/\code{base} are replaced by the swept grid values every cell exactly as the \code{k} argument replaces a supplied \code{leaf.prior}'s \code{k}; the object's other content - a \code{cgm} object's \code{split.probs}, a DART object's Dirichlet hyperparameters - rides every cell unchanged. Because \code{power}/\code{base} are grid axes here rather than ordinary scalars, they may be supplied alongside \code{tree.prior}; \code{dart} and \code{split.probs} would only duplicate what a supplied \code{tree.prior} already specifies, so combining either with it is an error naming both.
   }
}
\details{
  Crossvalidates \code{n.reps} replications against the crossproduct of given hyperparameter vectors \code{n.trees} \eqn{\times}{*} \code{k} \eqn{\times}{*} \code{power} \eqn{\times}{*} \code{base}. For each fit, either one fold is withheld as test data and \code{n.test - 1} folds are used as training data or \code{n * n.test} observations are withheld as test data and \code{n * (1 - n.test)} used as training. A replication corresponds to fitting all \eqn{K} folds in \code{"k-fold"} crossvalidation or a single fit with \code{"random subsample"}. The training data is used to fit a model and make predictions on the test data which are used together with the test data itself to evaluate the \code{loss} function.

  \code{loss} functions are either the default of average negative log-loss for binary outcomes and root-mean-squared error for continuous outcomes, misclassification rates for binary outcomes, or a \code{function} with arguments \code{y.test} and \code{y.test.hat}. \code{y.test.hat} is of dimensions equal to \code{length(y.test)} \eqn{\times}{*} \code{n.samples}, on the latent scale for binary outcomes. A third option is to pass a list of \code{list(function, evaluationEnvironment)}, so as to provide default bindings. The binary losses \code{"log"} and \code{"mcr"} cannot apply to continuous responses.

}
\value{
  An array with up to six dimensions, in order \code{rep} (length \code{n.reps}), \code{n.trees}, \code{k} (or \code{sd}), \code{power}, \code{base}, and \code{loss}. \code{rep} is always present. \code{n.trees}, \code{k}, \code{power}, and \code{base} are each omitted when \code{drop} is \code{TRUE} and the corresponding grid has length 1; with \code{drop = FALSE} they are always present, an absent \code{k} included, whose single cell is then named for the default it ran at. The trailing \code{loss} dimension, sized to however many values a single call to \code{loss} returns, is present only when that count is greater than one - independent of \code{drop} entirely; the default losses and an ordinary scalar-returning custom \code{loss} never contribute it. When none of the above survive, the result collapses to a plain vector of length \code{n.reps}. When the result remains an array, its \code{dimnames} name the swept values on each surviving hyperparameter axis - exact integer labels for \code{n.trees}, values rounded to 2 significant digits for the double-valued \code{k}/\code{power}/\code{base} axes, and the constructor call for a modelled \code{k} cell; the \code{loss} axis, when present, carries no per-slot names, and neither does \code{rep}.

  For method \code{"k-fold"}, each element is an average across the \eqn{K} fits. For \code{"random subsample"}, each element represents a single fit.

  The result is a bare array with no class, so the fit generics - \code{predict}, \code{extract}, \code{fitted}, \code{residuals} - do not apply to it; it is a table of losses, not a fit.
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

