\name{model.frame.dbartsData}
\alias{model.frame.dbartsData}
\alias{nobs.dbartsData}
\title{The Data of a dbartsData Object}
\description{
  Read the model's variables and the number of rows of a \code{\link{dbartsData}} object.
}
\usage{
\method{model.frame}{dbartsData}(formula, ...)
\method{nobs}{dbartsData}(object, ...)
}
\arguments{
  \item{formula, object}{A \code{\link{dbartsData}} object.}
  \item{\dots}{Not used.}
}
\details{
  \code{model.frame} returns the model's variables over the rows in, after \code{subset} and the \code{na.action}, as \code{model.frame} of an \code{lm} fit does: the response first, factors as factors, the weights in a \code{(weights)} column and the offsets in the \code{offset()} and \code{(offset)} columns, so that \code{\link{model.response}}, \code{\link{model.weights}} and \code{\link{model.offset}} read them. Its \code{"na.action"} attribute is the object's. On an object that keeps its frame it is made from the frame; where the frame holds a sparse column, which \code{model.frame} cannot read, it is the frame decoded from the block, as for an object with no frame.

  On an object with no frame (a sampler's data, \code{keep.data = FALSE}, the matrix interface, an environment) it is decoded from the block where that is exact: the response is returned in its recorded type and levels (a factor, a logical, a character vector, a \code{Surv} response where the object holds its time and status, the count matrix for a multinomial one) and named as the formula wrote it (\code{log(price)}, and \code{y} where no formula did), a factor predictor as a factor over its declared levels, and a matrix term such as \code{poly(x, 2)} as one matrix column again, named as the formula wrote it, with the column names the term gave. The block holds no more than the values, so a rejoined matrix has none of the attributes \code{poly()} or \code{scale()} give it (their class, coefficients and centers); every other column is as held.

  \code{nobs} counts every row in, rows of weight 0 included, which a BART fit holds in its trees; \code{nobs} of an \code{lm} fit leaves them out.
}
\value{
  A data frame, and a number.
}
\seealso{
  \code{\link{dbartsData}}
}
