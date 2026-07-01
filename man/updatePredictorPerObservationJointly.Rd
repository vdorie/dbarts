\name{updatePredictorPerObservationJointly}
\alias{updatePredictorPerObservationJointly}
\title{Jointly Update a Shared Predictor per Observation Across Samplers}
\description{
  Applies a \dQuote{partial} (per-observation) update of a single shared column to several
  \code{\linkS4class{dbartsSampler}} objects at once, installing an observation's new value only when
  it keeps every leaf node non-empty in every tree of every sampler.
}
\usage{
updatePredictorPerObservationJointly(samplers, x, column, updateState = NA)
}
\arguments{
  \item{samplers}{
    A \code{\linkS4class{dbartsSampler}}, or a list of \code{dbartsSampler} objects, that share the
    same (index-aligned) observations. The samplers are treated as peers; none is privileged.
  }
  \item{x}{
    A numeric vector of new values for the shared column, of length equal to the number of
    observations.
  }
  \item{column}{
    A single integer or character string identifying the shared column; it cannot be missing. The
    column is matched \emph{by name} across all samplers, so the shared variable may sit at different
    column positions in each. When given as an integer, it indexes the columns of the first sampler
    and the corresponding name is used to locate the column in the rest.
  }
  \item{updateState}{
    A logical determining if the local cache of each sampler's state should be updated after the
    update completes. If \code{NA}, each sampler's own \code{\link[=dbartsControl]{control}} default
    is used.
  }
}
\details{
  This is the multi-sampler analog of the \code{forceUpdate = "partial"} mode of
  \code{\link[=dbartsSampler-class]{setPredictor}}. A single sequential sweep installs each
  observation in every sampler at once, and only if its new value keeps every leaf non-empty in every
  tree of every chain of every sampler; otherwise that observation is rolled back to its previous
  value in all of them. Like the single-sampler partial mode, the update never changes tree structure
  -- it only re-routes observations.

  This supports embedding several BART components that condition on a common latent predictor (for
  example, a shared ability or trait in an item-response model) within a larger Gibbs/Metropolis
  sampler.
}
\value{
  A logical vector of length equal to the number of observations, \code{TRUE} where the observation's
  new value was installed and \code{FALSE} where it was rolled back to keep every tree valid. The
  same accept/reject decision is applied to every sampler, so an installed observation takes the new
  value in all of them; a rejected observation reverts to each sampler's own previous value. When the
  samplers agreed on the shared column before the call, they therefore continue to agree afterward.
}
\seealso{
  \code{\link{dbarts}}, \code{\linkS4class{dbartsSampler}}
}
