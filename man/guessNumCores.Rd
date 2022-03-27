\name{guessNumCores}
\title{Guess Number of Cores}
\alias{guessNumCores}
\description{
  Attempts to guess the number of CPU \sQuote{cores}, both physical and logical.
}
\usage{
guessNumCores(logical = FALSE)
}
\arguments{
  \item{logical}{
    A logical value. When \code{FALSE}, an estimate of the number of physical cores is returned. When \code{TRUE}, so-called \dQuote{logical} cores as also included.
  }
}
\details{
  Because of different definitions of cores used by different manufacturers, the distinction between logical and physical cores is not universally recognized. This function will attempt to use operating system definitions when available, which should usually match the CPU itself.
}
\value{
  An integer, or NA if no clear answer was obtained.
}
\author{
  Vincent Dorie: \email{vdorie@gmail.com}.
}
\keyword{parallel}
