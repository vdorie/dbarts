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
  \code{model.frame} returns the model's variables over the rows in, after \code{subset} and the \code{na.action}, as \code{model.frame} of an \code{lm} fit does: the response first, factors as factors, the weights in a \code{(weights)} column and the offsets in the \code{offset()} and \code{(offset)} columns, with a terms record, so that \code{\link{model.response}}, \code{\link{model.weights}} and \code{\link{model.offset}} read them. Its \code{"na.action"} attribute is the object's. On an object that keeps its frame it is made from the frame, a \code{poly()} term with its coefficients; a sparse column (a \code{sparseVector}, \code{\link{sparseFactor}} or \code{dgCMatrix}) comes back as stored, at the rows in, and is never made dense, so the frame holds classes a frame for \code{lm} could not.

  On an object with no frame (a sampler's data, \code{keep.data = FALSE}, the matrix interface, an environment) it is decoded from the block where that is exact: the response is returned in its recorded type and levels (a factor, a logical, a character vector, a \code{Surv} response where the object holds its time and status, the count matrix for a multinomial one) and named as the formula wrote it (\code{log(price)}, and \code{y} where no formula did), a factor predictor as a factor over its declared levels, a matrix term such as \code{poly(x, 2)} or \code{scale(x)} as one matrix column again, named as the formula wrote it, with the column names the term gave, and a column the block holds sparse as a sparse vector, sparse factor or sparse matrix. The columns come in the formula's order, those the block holds sparse after the others, and a column no term of the formula names (a hazard sampler's \code{period}) after them. The block holds no more than the values, so a rejoined matrix has none of the attributes \code{poly()} or \code{scale()} give it (their class, coefficients and centers); the terms record carries them in its \code{predvars}, as \code{lm}'s does. That record is the formula's, each sparse column listed as a term after the others, so that \code{lm} on the frame stops at it as it does on the frame a formula and data give; with no formula (the matrix interface), or where a hazard sampler added its \code{period}, the frame has none and \code{model.response} does not read it: its first column is the response.

  \code{nobs} counts the rows in with a positive weight, as the build counted them, evaluating nothing again, as \code{nobs} of a fit and of an \code{lm} fit do; rows of weight 0 stay in the object, and a BART fit holds them in its trees.
}
\value{
  A data frame, and a number.
}
\seealso{
  \code{\link{dbartsData}}
}
