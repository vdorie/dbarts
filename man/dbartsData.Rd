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
    missing = c("incorporate", "error"))
}
\arguments{
   \item{formula, data, test, subset, weights, offset, offset.test, factors, missing}{
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
    Factor level tables are retained so that test data are coded identically; a \code{sparseFactor} test column densifies to its integer codes over the training level table (implicit rows filled with the training reference level's code) before prediction, reusing the dense categorical test-code path.

    Data frame input is stored columnar: no \eqn{n \times p}{n x p} double matrix is retained. Code that previously reached into the data object's \code{x} slot expecting a plain matrix should call \code{\link{extract}(sampler, "predictors")} - or \code{as.matrix} on the slot - to obtain the numeric predictor matrix, factor columns as their integer codes.
  }
}
\value{
  An object of class \code{dbartsData}.
}
\seealso{
  \code{\link{dbarts}}, \code{\link{extract.dbartsSampler}}, \code{\link{sparseFactor}}
}
