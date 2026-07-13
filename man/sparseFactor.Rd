\name{sparseFactor}
\title{Sparse Unordered Factors}
\alias{sparseFactor}
\alias{sparseFactor-class}
\alias{show,sparseFactor-method}
\alias{length,sparseFactor-method}
\description{
  Represents an unordered factor sparsely: entries at the given positions carry the supplied levels, while every other position implicitly carries a designated \code{reference} level. Analogous to \code{Matrix::\link[Matrix]{sparseVector}}, whose implicit entry is zero.

  \strong{Not yet accepted by the samplers:} \code{\link{dbartsData}} recognizes a \code{sparseFactor} column but refuses to construct data from it with the error \dQuote{sparse categorical predictors are not yet supported}. The class exists so data sets can be assembled ahead of engine support; sparse \emph{ordinal} columns (\code{Matrix::sparseVector} or \code{dgCMatrix}) are unaffected and keep working.
}
\usage{
sparseFactor(x, levels, reference, i, length)
}
\arguments{
  \item{x}{
    The stored entries: a factor, a character vector, or integer level codes (which require an explicit \code{levels}). Missing values are not supported. When \code{i} is omitted, \code{x} is the complete (dense) vector and its non-\code{reference} entries become the stored ones.
  }
  \item{levels}{
    Character vector of factor levels. Defaults to \code{levels(x)} for a factor and \code{sort(unique(x))} for a character vector, the \code{\link{factor}} convention.
  }
  \item{reference}{
    The level every unstored position carries. Must be an element of \code{levels}; defaults to \code{levels[1]}, the baseline-contrast convention.
  }
  \item{i}{
    Optional integer vector of 1-based positions at which the entries of \code{x} sit, as in \code{Matrix::sparseVector}. When supplied, \code{x} and \code{i} must have equal length and \code{length} is required. Positions need not be sorted; duplicates are an error.
  }
  \item{length}{
    Total number of observations. Required alongside \code{i}; defaults to \code{length(x)} for dense input.
  }
}
\details{
  The constructor canonicalizes its input: entries are re-ordered by ascending position, entries whose level equals \code{reference} are dropped (they are the implicit value), and positions are stored 0-based in the \code{i} slot with 1-based level codes in the \code{values} slot. The remaining slots are \code{levels}, \code{reference}, and \code{length}.

  \code{show} prints the length, the number of stored entries, the level table, and the reference level. \code{length} returns the observation count, which lets a \code{sparseFactor} be a \code{\link{data.frame}} column.
}
\value{
  An object of S4 class \code{sparseFactor}.
}
\author{
  Vincent Dorie: \email{vdorie@gmail.com}.
}
\seealso{
  \code{\link{dbartsData}}
}
\examples{
f <- factor(c("a", "b", "a", "c", "a", "b"))
sf <- sparseFactor(f, reference = "a")
sf

## equivalently, supply only the non-reference entries and their positions
sparseFactor(c("b", "c", "b"), levels = c("a", "b", "c"), reference = "a",
             i = c(2L, 4L, 6L), length = 6L)
}
\keyword{classes}
