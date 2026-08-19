\name{dbartsMixedMatrix}
\title{Mixed Dense/Sparse Predictor Container}
\alias{dbartsMixedMatrix}
\alias{dim.dbartsMixedMatrix}
\alias{dimnames.dbartsMixedMatrix}
\alias{[.dbartsMixedMatrix}
\alias{as.matrix.dbartsMixedMatrix}
\description{
  The container a data frame's predictor columns assemble into whenever any
  column is sparse-backed (a \code{\link{sparseFactor}} column, or a
  \code{Matrix::sparseVector}/\code{dgCMatrix} column assigned into the
  frame): dense columns stay a per-column list and sparse columns stay one
  shared \code{dgCMatrix}, rather than materializing one \eqn{n \times
  p}{n x p} dense matrix. \code{\link{makeModelMatrixFromDataFrame}} returns
  one whenever \code{x} carries a sparse column; \code{\link{dbartsData}}'s
  x/y interface (and so \code{\link{dbarts}}'s and \code{\link{bart2}}'s)
  returns one under the same condition and stores it in the resulting
  \code{dbartsData} object's \code{x}/\code{x.test} slot. The formula
  interface refuses a sparse column outright (a bare S4 column does not
  survive \code{\link{model.frame}}), so a \code{dbartsMixedMatrix} only
  ever arises from one of those two entrances.
}
\usage{
\method{dim}{dbartsMixedMatrix}(x)
\method{dimnames}{dbartsMixedMatrix}(x)
\method{[}{dbartsMixedMatrix}(x, i, j, drop = TRUE)
\method{as.matrix}{dbartsMixedMatrix}(x, \dots)
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

  What is \strong{not} supported: matrix-style subassignment (\code{x[i, j]
  <- value} errors), and coercion to a \code{\link{data.frame}}
  (\code{as.data.frame(x)} errors - there is no \code{as.data.frame}
  method). \code{is.matrix(x)} is \code{FALSE}. Code that needs either
  should call \code{as.matrix(x)} first.
}
\value{
  \code{dim} returns an integer vector of length two. \code{dimnames}
  returns a length-two list, the first element always \code{NULL}. \code{[}
  with \code{j} supplied returns a numeric vector or matrix, matching
  \code{drop}; with \code{j} missing it returns another
  \code{dbartsMixedMatrix}. \code{as.matrix} returns a plain numeric
  matrix.
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
