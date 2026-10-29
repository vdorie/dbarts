\name{dbartsMixedMatrix}
\title{Mixed Dense/Sparse Predictor Container}
\alias{dbartsMixedMatrix}
\alias{dim.dbartsMixedMatrix}
\alias{dimnames.dbartsMixedMatrix}
\alias{[.dbartsMixedMatrix}
\alias{as.matrix.dbartsMixedMatrix}
\alias{as.data.frame.dbartsMixedMatrix}
\description{
  The container a data frame's predictor columns assemble into whenever any
  column is sparse-backed (a \code{\link{sparseFactor}} column, or a
  \code{Matrix::sparseVector}/\code{dgCMatrix} column assigned into the
  frame): dense columns stay a per-column list and sparse columns stay one
  shared \code{dgCMatrix}, rather than materializing one \eqn{n \times
  p}{n x p} dense matrix. \code{\link{makeModelMatrixFromDataFrame}} returns
  one whenever \code{x} carries a sparse column; \code{\link{dbartsData}}'s
  x/y interface (and so \code{\link{dbarts}}'s and \code{\link{bart}}'s)
  returns one under the same condition and stores it in the resulting
  \code{dbartsData} object's \code{x}/\code{x.test} slot. A fit or sampler
  made from a formula or a data frame holds one as its \code{data@x}
  whether or not a column is sparse. The formula
  interface accepts a sparse column too, named like any other (including
  through \code{.}): a bare S4 column cannot survive
  \code{\link{model.frame}}, so it is lifted out of \code{data} first and
  re-attached, row-subset under \code{subset} and \code{na.action}, to the
  assembled predictor matrix afterward.
}
\usage{
\method{dim}{dbartsMixedMatrix}(x)
\method{dimnames}{dbartsMixedMatrix}(x)
\method{[}{dbartsMixedMatrix}(x, i, j, drop = TRUE)
\method{as.matrix}{dbartsMixedMatrix}(x, \dots)
\method{as.data.frame}{dbartsMixedMatrix}(x, row.names = NULL, optional = FALSE, \dots)
}
\arguments{
  \item{x}{
    An object of class \code{dbartsMixedMatrix}.
  }
  \item{i, j}{
    Row and column indices, as for a matrix. Supplying \code{j} densifies
    (see \sQuote{Details}); \code{i} alone keeps the container.
  }
  \item{drop}{
    Logical; as for a matrix, only consulted when \code{j} is supplied.
  }
  \item{row.names}{
    \code{NULL}, the default, for a data frame's default row names, or a
    vector of row names.
  }
  \item{optional}{
    Unused.
  }
  \item{\dots}{
    Unused.
  }
}
\details{
  \code{dim} and \code{dimnames} - and so the base \code{nrow}, \code{ncol},
  and \code{colnames}, which are defined in terms of them and need no
  method of their own - read like a matrix's; row names are always
  \code{NULL}. \code{x[, j]} (or \code{x[, "name"]}) returns the
  indexed/named column(s) as an ordinary numeric vector or matrix.
  \code{x[i, ]} keeps the container, row-subset. \code{x[i]} with \code{j}
  omitted also row-subsets, unlike a plain matrix's single-index linear
  indexing. \code{as.matrix} materializes the whole object to a dense
  \eqn{n \times p}{n x p} numeric matrix (a factor column's codes,
  zero-based) and is the escape hatch for any operation the container
  itself does not support.

  \code{as.data.frame} gives the predictors as the model sees them, one
  column per column of the container, named by \code{colnames(x)}: numeric
  columns as numeric, and a factor column as a factor over its declared
  levels, which keeps a level no row takes. For a fit made from a formula of
  plain variables, \code{predict} and \code{setTestPredictor} take that data
  frame as new data. A formula that transforms or expands a variable
  (\code{I(z^2)}, \code{factor(i)}, a matrix term, \code{poly(z, 2)}) names
  columns that are not the variables it reads, so new data for it is a data
  frame of those variables, not this one.

  What is \strong{not} supported: matrix-style subassignment (\code{x[i, j]
  <- value} errors). \code{is.matrix(x)} is \code{FALSE}, and a base function
  that does not dispatch on the class, such as \code{cbind}, does not see a
  matrix. Code that needs a matrix should call \code{as.matrix(x)} first.
}
\value{
  \code{dim} returns an integer vector of length two. \code{dimnames}
  returns a length-two list, the first element always \code{NULL}. \code{[}
  with \code{j} supplied returns a numeric vector or matrix, matching
  \code{drop}; with \code{j} missing it returns another
  \code{dbartsMixedMatrix}. \code{as.matrix} returns a plain numeric
  matrix, and \code{as.data.frame} a data frame.
}
\seealso{
  \code{\link{makeModelMatrixFromDataFrame}}, \code{\link{dbartsData}},
  \code{\link{sparseFactor}}
}
\examples{
df <- data.frame(a = 1:5)
df$b <- Matrix::sparseVector(c(1, 2), c(1, 3), length = 5)
mm <- makeModelMatrixFromDataFrame(df)

class(mm)
dim(mm)
colnames(mm)
mm[, "a"]
as.matrix(mm)
}
\keyword{classes}
