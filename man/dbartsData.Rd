\name{dbartsData}
\alias{dbartsData}
\title{Discrete Bayesian Additive Regression Trees Sampler Data}
\description{
  Convenience function to create a data object for use with a \code{\link{dbarts}} sampler.
}
\usage{
dbartsData(
    formula, data, test, subset, weights,
    offset, offset.test = offset,
    factors = c("categorical", "indicators"),
    na.action = dbarts::na.keepPredictors,
    bases = NULL,
    counts = NULL)
}
\arguments{
   \item{formula, data, test, subset, weights, offset, offset.test, factors, na.action}{
     As in \code{\link{dbarts}}. Retains backwards compatibility with \code{\link{bart}}, so that \code{formula}/\code{data} can be a \code{\link{formula}}/\code{\link{data.frame}} pair, or a pair of \code{x.train}/\code{y.train} matrices/vector. \code{offset}/\code{offset.test} take the matrix interface's shape from whether \code{counts} is also supplied: with \code{counts}, an \code{n} x \code{K} (\code{nTest} x \code{K}) numeric matrix installs the per-category shift described under \code{counts} below, and a flat vector is refused; without \code{counts}, only a flat vector (or single number) is accepted, and a matrix is refused. \code{formula} may itself be a \code{dbartsData} object, in which case it is returned as-is; any other argument supplied alongside it is then ignored, with a warning (class \code{dbartsIgnoredArgWarning}).
   }
   \item{bases}{
     Optional list of the per-forest amplitude bases a multi-forest model combines its forests through, one element per forest and \code{NULL} where a forest carries the implicit intercept its single amplitude scales. Each non-null element is a numeric matrix with one row per observation and one amplitude per column, restricted by \code{subset} exactly as \code{weights} is. This is the data half of a multi-forest specification, and is normally written by \code{\link{dbarts}}'s \code{forests = } argument rather than here; supplying it directly selects the multi-forest model with the engine's default per-forest configuration, and \code{\link{forest}} declarations may then configure the forests one at a time.
   }
   \item{counts}{
     Optional \code{n} x \code{K} matrix of non-negative integer counts selecting the multinomial (softmax) model, column \code{k} holding category \code{k}'s successes and every row carrying at least one trial. It is restricted by \code{subset} exactly as \code{weights} is, and it IS the response: \code{y} is derived as its row sums - the trials \eqn{n_i}{n_i} - so no separate response need be supplied, and any that is supplied is replaced. At least two categories are required. Column names, when present, label the categories. \code{\link{dbarts}(family = "multinomial")} writes this argument from a count-matrix or factor response rather than requiring it here. Paired with \code{counts}, a matrix \code{offset}/\code{offset.test} is the per-category shift added to the raw per-category values BEFORE the softmax - which is the only place a shift means anything under one, since the blend is invariant to a common per-observation shift and a flat shift would be added after it. Only the row-centred part is identified: adding a constant to a whole row leaves every reported probability unchanged.
   }
}
\details{
  \subsection{Data frame ingestion}{
    Data frame columns - through the formula interface or a frame passed as \code{x.train} - map to predictors as follows:
    \itemize{
      \item numeric and logical columns enter as ordinal predictors;
      \item unordered factors each become a single categorical predictor under the default \code{factors = "categorical"}, with splits sending a subset of levels down each branch; \code{factors = "indicators"} instead expands each factor into binary indicator columns, as previous versions always did - a factor past 100 levels builds this expansion as a sparse (\code{dgCMatrix}-backed) block automatically, with no user-facing argument: the tree structure grown is the same as the dense expansion's, though leaf-mean draws can differ from it in the last bit;
      \item ordered factors enter as their own predictor kind, coded numerically by level order and split at a threshold on those codes; the split grid is the midpoints between consecutive \emph{declared} levels, one candidate per level boundary, so \code{n.cuts} does not apply to such a column and adjacent levels are always separable;
      \item matrix columns (e.g. from \code{\link{poly}}) splice in as ordinal predictors, one per column;
      \item \code{Matrix::sparseVector} and \code{dgCMatrix} columns (assigned into the frame; they do not survive \code{data.frame(...)} or \code{I()}) enter as sparse ordinal predictors. Both interfaces accept them for training: a bare S4 column does not survive \code{\link{model.frame}}, so the formula path lifts any such column out of \code{data} before \code{model.frame} runs and re-attaches it, row-subset under \code{subset} and \code{na.action} alike, to the assembled predictor matrix afterward - named in \code{formula} like any other column, including through \code{.}, with no marker. A sparse name wrapped in \code{poly()}, \code{ns()}, \code{log()}, \code{offset()}, or a \code{':'}/\code{'*'} interaction is refused by name, as an interaction among dense terms already is. \code{test}/\code{x.test} accepts them under either interface, since a test set is never run through \code{model.frame};
      \item \code{\link{sparseFactor}} columns enter as a single sparse categorical predictor, coded and binned identically to a dense factor of the same values. Both interfaces accept them, the same pull-out/re-attach the other sparse column kinds get on the formula path.
    }
    On the formula path, a sparse column re-attaches after every dense term regardless of where it was named, each in the order its name appears in \code{formula} (after \code{.} expansion) rather than its position in \code{data}; \code{test} is matched to \code{x.train} by name (or by position when neither carries names), so the two stay aligned regardless of this reordering.
    Factor level tables are retained so that test data are coded identically, and a factor predictor's level count - ordered or unordered - is the \emph{declared} level table rather than the levels the training rows happen to take. A level declared but unobserved in training therefore keeps its own bin, and test or replacement data carrying it are accepted at every entrance (creation, \code{setTestPredictor}, \code{predict}, \code{setPredictor}); conversely a value that is not an existing integer level code is refused at those same entrances, naming the kind, for ordered and unordered factors alike. An unordered factor is limited to 65535 levels and an ordered one to 65534 - one lower, since an ordered factor spends a code on the upper bin of its \eqn{K - 1}{K - 1} cut grid - and a wider one is an error naming the column. This holds equally for dense factors, \code{sparseFactor} columns, and dense factors inside a mixed dense/sparse frame. Drop unused levels before fitting (\code{\link{droplevels}}) if you would rather they not be modeled. \code{test} accepts the same column types as \code{x.train}, including \code{sparseFactor} and sparse ordinal columns. A \code{sparseFactor} or sparse ordinal test column recodes over the training level table and stays resident - through creation and \code{setTestPredictor} - rather than densifying at ingestion, the same storage tier rule as training. \code{predict} and \code{getTrees(newdata = )} code a test set the same way and then route its rows through the trees off that storage, materializing no dense matrix of their own.

    Data frame input is stored columnar: no \eqn{n \times p}{n x p} double matrix is retained. Code that previously reached into the data object's \code{x} slot expecting a plain matrix should call \code{\link{extract}(sampler, "predictors")} - or \code{as.matrix} on the slot - to obtain the numeric predictor matrix, factor columns as their integer codes.
  }
}
\value{
  An object of class \code{dbartsData}.
}
\seealso{
  \code{\link{dbarts}}, \code{\link{extract.dbartsSampler}}, \code{\link{sparseFactor}}
}
\examples{
set.seed(3)
n <- 60L
df <- data.frame(
  x = runif(n),
  g = factor(sample(c("a", "b", "c"), n, replace = TRUE))
)
df$y <- df$x + 0.5 * (df$g == "b") + rnorm(n, 0, 0.1)

## the formula interface; 'g' enters as one categorical predictor
data <- dbartsData(y ~ x + g, df)
data

## the same data through the backwards-compatible x/y interface, with a
## test set that is coded against the training level table
train <- df[1:40, ]
test  <- df[41:60, c("x", "g")]
data <- dbartsData(train[, c("x", "g")], train$y, test)

## build a sampler on it, rather than letting dbarts() ingest again
sampler <- dbarts(data, control = dbartsControl(
  n.chains = 1L, n.threads = 1L, n.samples = 10L, n.burn = 10L,
  n.trees = 25L, updateState = FALSE))
samples <- sampler$run()
str(samples$test)
}
