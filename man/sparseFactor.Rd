\name{sparseFactor}
\title{Sparse Unordered Factors}
\alias{sparseFactor}
\alias{sparseFactor-class}
\alias{show,sparseFactor-method}
\alias{length,sparseFactor-method}
\alias{[,sparseFactor-method}
\alias{as.character,sparseFactor-method}
\alias{format.sparseFactor}
\alias{as.data.frame.sparseFactor}
\alias{[[,sparseFactor-method} % ]]
\alias{[<-,sparseFactor-method}
\alias{c,sparseFactor-method}
\alias{as.integer,sparseFactor-method}
\alias{[[<-,sparseFactor-method} % ]]
\alias{levels<-,sparseFactor-method}
\alias{rep,sparseFactor-method}
\alias{droplevels.sparseFactor}
\alias{as.vector,sparseFactor-method}
\alias{is.na,sparseFactor-method}
\alias{xtfrm,sparseFactor-method}
\alias{length<-,sparseFactor-method}
\alias{Ops,sparseFactor,ANY-method}
\alias{Ops,ANY,sparseFactor-method}
\alias{Ops,sparseFactor,sparseFactor-method}
\alias{str.sparseFactor}
\alias{summary.sparseFactor}
\alias{unique.sparseFactor}
\alias{duplicated.sparseFactor}
\description{
  Represents an unordered factor sparsely: entries at the given positions carry the supplied levels, while every other position implicitly carries a designated \code{reference} level. Analogous to \code{Matrix::\link[Matrix]{sparseVector}}, whose implicit entry is zero.

  Accepted as a predictor by \code{\link{dbartsData}} through either interface: a data frame passed as \code{x.train} may mix dense numeric/factor columns with sparse ordinal columns (\code{Matrix::sparseVector} or \code{dgCMatrix}) and \code{sparseFactor} columns in any combination, and a \code{formula}'s \code{data} accepts one too, named like any other column (including through \code{.}): a bare S4 column cannot survive \code{\link{model.frame}}, so \code{dbartsData} lifts it out of \code{data} first and re-attaches it, row-subset under \code{subset} and \code{na.action}, to the assembled predictor matrix afterward. A \code{sparseFactor} column enters as one categorical predictor and bins bitwise-identically to a dense factor of the same values.

  A test frame (\code{test}, or \code{newdata} to \code{predict}) is read the way a dense-factor frame is: the model's terms are replayed on its dense columns (so a transformed term such as \code{log(z)} works and extra columns, the response included, are dropped) and its \code{sparseFactor} columns are used by name. The same holds for a test set (\code{test}/\code{x.test}): a \code{sparseFactor} test column is recoded over the training level table and stays resident - through creation and \code{setTestPredictor} - rather than densifying at ingestion. \code{predict} and \code{getTrees(newdata = )} code it the same way and then route its rows through the trees off that storage, materializing no dense matrix of their own.

  Mixing a \code{sparseFactor} (or any sparse ordinal) column with ordinary dense columns changes how a default starting \code{sigma} is estimated. The unmixed case fits an \code{\link{lm}} on the training design and uses its residual standard deviation; as soon as \emph{any} column is sparse-backed, the design is not densified for this purpose and the estimate falls back to the unconditional \code{sd(y)} instead - discarding every dense column's information along with the sparse one's. The fallback warns (class \code{dbartsSparseSigmaFallbackWarning}) rather than passing silently. This is current behavior, not a documented guarantee; supply \code{\link{dbarts}}'s \code{sigma} (\code{sigest} in \code{\link{bart}}) explicitly if the default matters to you.
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
    The level every unstored position carries. Must be an element of \code{levels}; defaults to \code{levels[1]}, the baseline-contrast convention. Storage holds only the non-reference entries, so choosing the \strong{most common level} as the reference maximizes the memory win; codes and every draw are identical under any choice of reference, so this is a storage tuning knob, not a modeling one.
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

  A \code{sparseFactor} answers the operations a data frame reaches a factor column with as a \code{\link{factor}} does. \code{x[i]} subsets by position (positive, negative, zero, logical and repeated indices) and returns a \code{sparseFactor} over the same levels, mapping the stored positions without densifying; \code{drop = TRUE} drops the levels no selected row takes, and with no row selected keeps the reference level alone, since a \code{sparseFactor} needs one level. \code{x[i] <- value} and the double-bracket assignment take labels that are levels of \code{x} (as a character vector, factor or \code{sparseFactor}) and may extend the vector when no position is left unassigned; \code{levels(x) <- value} renames the levels, the reference with them, merging repeated names as it does for a factor. \code{c(x, ...)} combines with factors and \code{sparseFactor}s over the union of their levels. \code{as.character}, \code{as.vector} and \code{format} give the level labels, \code{as.integer} the level codes, and \code{levels} (read from the class's slot), \code{is.na} (always \code{FALSE}), \code{xtfrm} (so \code{order} and \code{sort}), \code{unique}, \code{duplicated}, \code{rep}, \code{droplevels}, \code{summary}, \code{==} and \code{!=} (which refuse two factors with different level sets, as for factors) behave as for a factor, and \code{str} prints a factor's line. \code{factor}, \code{as.factor} and \code{table} read it through those methods, over the levels present. Together these let a data frame holding one be subset by row, assigned into, printed and \code{str}-ed. \code{rbind} of such frames returns the column as an ordinary factor; it needs R 4.6.0 or later, since earlier versions of base R's \code{match} refuse an S4 object (convert the column with \code{factor} first).

  A \code{sparseFactor} cannot hold a missing value, so where a factor would return \code{NA} it is refused: an \code{NA} or out-of-range row index, a replacement value that is not a level, an extension past a gap, a level named \code{NA}, and \code{length<-}. A character index has no names to match and is refused. Not supported: \code{complete.cases} on a data frame holding one (use \code{\link{na.omit}}, or \code{is.na} per column), and \code{c(f, x)} with a factor \code{f} first, which \code{c.factor} answers with a list (put the \code{sparseFactor} first). Also not supported: \code{relevel}, \code{as.numeric} (\code{as.integer} gives the level codes), \code{rep_len} and \code{rep.int} (\code{rep} works), and the \code{exclude} argument of \code{droplevels}, which is ignored; \code{droplevels} on a data frame leaves a \code{sparseFactor} column untouched, as base R touches only factor columns, so call it on the column.

  A \code{Matrix::sparseVector} or \code{dgCMatrix} column subsets in a data frame through Matrix's own method, but printing such a frame shows a placeholder for the column and warns of a corrupt data frame; that formatting belongs to Matrix. Use \code{as.numeric} on the column to display it.
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
