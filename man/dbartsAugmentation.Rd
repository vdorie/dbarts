\name{dbartsAugmentation}
\alias{dbartsAugmentation}
\alias{dbartsDrawLatents}
\alias{dbartsWorkingResponse}
\title{Augmentation Draws For A Composed Sampler}
\description{
  The per-observation data augmentation a non-gaussian \code{\linkS4class{dbartsSampler}} runs
  inside each sweep, callable from R against R's own random number stream.
  \code{dbartsDrawLatents} draws the augmentation variable at a supplied fit - the quantity
  \code{\link[=dbartsSampler-class]{getLatents}} reports for that family - and
  \code{dbartsWorkingResponse} turns a drawn latent into the response a host regresses on. Neither
  touches a sampler: a host writing its own Gibbs step around BART, or around anything else, need
  not reimplement the family's augmentation.
}
\usage{
dbartsDrawLatents(family, fit, y, weights = NULL, offset = NULL, sigma = NULL,
                  dispersion = NULL, cutpoints = NULL, df = NULL)

dbartsWorkingResponse(family, latent, y, weights = NULL, offset = NULL,
                      dispersion = NULL)
}
\arguments{
  \item{family}{
    One of \code{"probit"}, \code{"logistic"}, \code{"ordinal"}, \code{"aft"}, \code{"nbinom"} or
    \code{"student"}, matched by \code{\link{match.arg}}; there is no default. \code{"student"} is
    the Student-t residual distribution (\code{resid.dist = student()} in \code{\link{dbarts}}), a
    gaussian response under a scale mixture.
  }
  \item{fit}{
    A numeric vector of each observation's location \emph{without} the offset -
    \eqn{f(x_i)}{f(x_i)}, what \code{\link[=dbartsSampler-class]{getFitsWithoutOffset}} reports, or
    however the host computes it. Its length sets \eqn{n}{n}.
  }
  \item{latent}{
    A numeric vector of length \eqn{n}{n} of drawn latents, as \code{dbartsDrawLatents} returns;
    positive for \code{"logistic"} and \code{"nbinom"}, whose latents are precisions.
  }
  \item{y}{
    The response, of length \eqn{n}{n}, constrained exactly as the sampler constrains it: 0 or 1 for
    \code{"probit"} and \code{"logistic"}, an integer category index in \eqn{[1, K]}{[1, K]} for
    \code{"ordinal"}, a non-negative integer count for \code{"nbinom"}. For \code{"aft"} it is the
    log observed (censoring) time, the bound the drawn log survival time is truncated below at;
    \code{"student"} constrains nothing.
  }
  \item{weights}{
    \code{"logistic"} only, refused by name elsewhere: observation counts, positive whole numbers,
    each entering as the sum of that many independent \eqn{PG(1, \psi)}{PG(1, psi)} draws.
  }
  \item{offset}{
    A numeric vector of length \eqn{n}{n} added to \code{fit} to form the linear predictor, or
    \code{NULL} for none.
  }
  \item{sigma}{
    \code{"aft"} and \code{"student"} only: the residual standard deviation, in response units.
  }
  \item{dispersion}{
    \code{"nbinom"} only, and required there: the current dispersion \eqn{r}{r}, a positive whole
    number, as \code{\link[=dbartsSampler-class]{getDispersion}} reports it.
  }
  \item{cutpoints}{
    \code{"ordinal"} only, and required there: the \eqn{K - 1}{K - 1} strictly increasing category
    boundaries in force.
  }
  \item{df}{
    \code{"student"} only, and required there: the degrees of freedom \eqn{\nu}{nu}.
  }
}
\details{
  \strong{The offset convention.} \code{fit} is the location \emph{without} the offset, the engine's
  own convention at every augmentation site, and the linear predictor is formed internally as
  \code{fit + offset}. The two arguments are therefore not interchangeable, and the pair of
  functions is what tells them apart: \code{dbartsDrawLatents(fit = f, offset = o)} and
  \code{dbartsDrawLatents(fit = f + o, offset = NULL)} draw from the same law and agree, while
  \code{dbartsWorkingResponse} \emph{differs} between the two, the working response being stated in
  units the offset has been taken out of. Pass the offset as \code{offset} to both.

  \strong{Why two functions.} For \code{"probit"}, \code{"ordinal"} and \code{"aft"} the latent is a
  LOCATION and the working response is that location less the offset. For \code{"logistic"} and
  \code{"nbinom"} it is a PRECISION - the Polya-Gamma \eqn{\omega_i}{omega_i} - and what a host
  regresses on is \eqn{w_i (y_i - 1/2) / \omega_i - o_i}{w_i (y_i - 1/2) / omega_i - o_i} and
  \eqn{(y_i - r) / (2 \omega_i) - o_i}{(y_i - r) / (2 omega_i) - o_i} respectively, under per-row
  precision weights \eqn{\omega_i}{omega_i}. For \code{"student"} the latent is the scale-mixing
  \eqn{\lambda_i}{lambda_i}, which weights the row rather than entering the response, so the working
  response is the gaussian one. A host that regresses on the latent itself is fitting a different
  model, which is the composition defect these functions exist to prevent.

  \strong{Random numbers.} The draws come from R's own stream, through a generator created and
  destroyed within the call: \code{\link{set.seed}} reproduces them, they advance the stream, and
  they compose with any other R draw. They never touch a sampler's chain generator, which is seeded
  separately and never advances R's stream. That stream is \emph{not} thread safe, so these
  functions must not be called from a worker thread; call them on the main thread and hand the
  result to the workers.
}
\value{
  A numeric vector of length \eqn{n}{n}. \code{dbartsDrawLatents} carries a \code{"quantity"}
  attribute of \code{"location"} (\code{"probit"}, \code{"ordinal"}, \code{"aft"}) or
  \code{"precision"} (\code{"logistic"}, \code{"nbinom"}, \code{"student"}), matching what
  \code{\link[=dbartsSampler-class]{getLatents}} reports for that family.
}
\seealso{
  \code{\link[=dbartsSampler-class]{dbartsSampler}}, whose \code{getLatents},
  \code{getFitsWithoutOffset} and \code{getDispersion} methods report the same quantities from
  inside a running sampler.
}
\examples{
set.seed(99)
f <- rnorm(50L)
o <- rep(0.25, 50L)
y <- rbinom(50L, 1L, pnorm(f + o))

# one Albert-Chib step: the truncated-normal latent at the current fit, then
# the response a host regresses on
z <- dbartsDrawLatents("probit", fit = f, y = y, offset = o)
attr(z, "quantity")
head(dbartsWorkingResponse("probit", latent = z, y = y, offset = o))

# logistic draws a PRECISION, so its working response is a quotient, not the
# latent shifted
w <- dbartsDrawLatents("logistic", fit = f, y = y, offset = o)
head(dbartsWorkingResponse("logistic", latent = w, y = y, offset = o))
}
\keyword{models}
