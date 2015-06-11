\name{guessNumCores}
\title{Guess Number of Cores}
\alias{guessNumCores}
\alias{makeind}
\description{
  Attempts to guess the number of CPU \sQuote{cores}, both physical and logical.
}
\usage{
  guessNumCores(logical = false)
}
\arguments{
  \item{logical}{Logical. When \code{FALSE}, estimate of the number of physical cores is returned.}
}
\details{
  Because of different definitions of cores used by different manufacturers, the
  distinction between logical and physical cores is not universally recognized. This
  function will attempt to use operating system definitions when available, which should
  usually match the CPU itself.
}
\value{
  An integer, or NA if no clear answer was obtained.
}
\author{
  Vincent Dorie: \email{vjd4@nyu.edu}.
}
\keyword{parallel}
