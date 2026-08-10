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
    missing = c("incorporate", "error"),
    treatment = NULL)
}
\arguments{
   \item{formula, data, test, subset, weights, offset, offset.test, factors, missing, treatment}{
     As in \code{\link{dbarts}}. Retains backwards compatibility with \code{\link{bart}}, so that \code{formula}/\code{data} can be a \code{\link{formula}}/\code{\link{data.frame}} pair, or a pair of \code{x.train}/\code{y.train} matrices/vector.
   }
}
\details{
  \subsection{Data frame ingestion}{
    Data frame columns - through the formula interface or a frame passed as \code{x.train} - map to predictors as follows:
    \itemize{
      \item numeric and logical columns enter as ordinal predictors;
      \item unordered factors each become a single categorical predictor under the default \code{factors = "categorical"}, with splits sending a subset of levels down each branch; \code{factors = "indicators"} instead expands each factor into binary indicator columns, as previous versions always did;
      \item ordered factors enter as ordinal predictors, coded numerically by level order;
      \item matrix columns (e.g. from \code{\link{poly}}) splice in as ordinal predictors, one per column;
      \item \code{Matrix::sparseVector} and \code{dgCMatrix} columns (assigned into the frame; they do not survive \code{data.frame(...)} or \code{I()}) enter as sparse ordinal predictors;
      \item \code{\link{sparseFactor}} columns enter as a single sparse categorical predictor, coded and binned identically to a dense factor of the same values. Only the x/y interface accepts them (a bare S4 column cannot survive \code{\link{model.frame}}, so the formula path refuses one explicitly).
    }
    Factor level tables are retained so that test data are coded identically, and a categorical predictor's category count is the \emph{declared} level table rather than the levels the training rows happen to take. A level declared but unobserved in training therefore keeps its own bin, and test or replacement data carrying it are accepted at every entrance (creation, \code{setTestPredictor}, \code{predict}, \code{setPredictor}); this holds equally for dense factors, \code{sparseFactor} columns, and dense factors inside a mixed dense/sparse frame. Drop unused levels before fitting (\code{\link{droplevels}}) if you would rather they not be modeled. \code{test} accepts the same column types as \code{x.train}, including \code{sparseFactor} and sparse ordinal columns. A \code{sparseFactor} or sparse ordinal test column recodes over the training level table and stays resident - through creation and \code{setTestPredictor} - rather than densifying at ingestion, the same storage tier rule as training. \code{predict} and \code{getTrees(newdata = )} code a test set the same way and then materialize it to a dense matrix at the call boundary, since those entry points evaluate already-built trees on resolved values.

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
