\name{forest}
\alias{forest}
\title{Forest Specification for Multi-Forest Models}
\description{
  Build the specification of one forest of a multi-forest model. Pass a list of
  these as the \code{forests} argument of \code{\link{dbarts}} or
  \code{\link{dbartsSpec}}, which fits the mean as a weighted sum of the
  declared ensembles. Every knob is per forest and carried on the one
  constructor, so the fitting functions grow exactly one argument however many
  forests a model has.
}
\usage{
forest(
    basis = NULL, vars = NULL,
    n.trees = NULL, base = NULL, power = NULL, sd = NULL,
    interactions = NULL, blocks = NULL,
    amplitude.prior.variance = NULL, update.amplitude = NULL)
}
\arguments{
  \item{basis}{
    The data this forest's amplitudes multiply: a one-sided formula
    (\code{~ factor(z)}), evaluated against the fit's own \code{data} and then
    in its own environment, or an already-evaluated vector or matrix. It
    expands by R's model-matrix rule - a factor becomes its level indicators,
    one amplitude per level, with no reference level dropped, since the forest
    carries no intercept of its own - so the forest enters the mean as
    \eqn{(\sum_k a_k B_k(x_i)) f(x_i)}{(sum_k a_k B_k(x_i)) f(x_i)}. Any
    forest may carry one, of any width: a two-level factor gives the pair whose
    amplitudes are \eqn{(b_0, b_1)}{(b0, b1)}, a wider factor one amplitude per
    level, and a numeric vector or matrix is already those columns. Level order
    is meaningful: the second level is the one \eqn{b_1}{b1} scales. A forest
    that declares none takes the implicit intercept its single amplitude
    \eqn{a}{a} scales; every forest past the first needs one, since the
    amplitudes multiplying it are what distinguish it from the first.
  }
  \item{vars}{
    Optional restriction of this forest to a subset of the model matrix, by
    column name or index; \code{NULL} leaves it reading every predictor. Any
    forest may be restricted.
  }
  \item{n.trees, base, power}{
    This forest's tree count and tree-structure prior. \code{NULL} takes the
    engine's default, which for a basis forest is 50 trees at
    \code{base = 0.25}, \code{power = 3} - shallower and fewer than a
    prognostic forest's, the modulating surface normally being the smoother of
    the two. On the FIRST forest these are the fit's own \code{control@n.trees}
    and \code{tree.prior}, which they restate rather than add to. When a value
    here disagrees with an explicitly supplied \code{control@n.trees} or
    tree-prior \code{base}/\code{power}, this one governs the fit, being the
    more specific of the two declarations.
  }
  \item{sd}{
    This forest's prior scale, in units of the response family's own latent
    scale and per unit of basis row norm: the total \eqn{a\,f(x)}{a f(x)} or
    \eqn{(b_1 - b_0)f(x)}{(b1 - b0) f(x)} is placed at \code{sd} of them. The
    unit is \eqn{\mathrm{sd}(y)}{sd(y)} for a gaussian response, \code{1} for
    \code{"probit"} and \eqn{\pi/\sqrt{3}}{pi/sqrt(3)} for \code{"logistic"} -
    the standard deviation of the link's own error law, a latent model having
    no response standard deviation to name. A forest whose \code{basis} rows
    have median non-zero norm \eqn{c}{c} contributes the scale named here, the
    calibration map dividing \eqn{c}{c} out, so rescaling a basis column does
    not silently rescale the prior. Which of the two channels carries it depends
    on whether the forest has a \code{basis}: without one it is the half-Cauchy
    median of the forest's scalar amplitude, with one it scales the forest's own
    node prior. The two channels default differently, and neither default is a
    bare constant. A forest with NO basis takes \code{2} under a gaussian
    response, where the unit is the response's own
    \eqn{\mathrm{sd}(y)}{sd(y)} and a drawn \eqn{\sigma}{sigma} absorbs the
    difference, and \code{1} under \code{"probit"} and \code{"logistic"}, where
    the unit is the link's fixed error scale and nothing does. A forest WITH a
    basis takes \eqn{\sqrt{2/K}}{sqrt(2/K)} in a model of \eqn{K}{K} forests, so
    that declaring more of them does not widen the prior on the combined
    location without bound. \eqn{K = 2}{K = 2} is the fixed point of both
    statements, \eqn{\sqrt{2/2} = 1}{sqrt(2/2) = 1}. A value declared here
    overrides its default and keeps its per-forest reading at every \eqn{K}{K},
    so \code{sd = 1} on each basis forest recovers the pre-\eqn{K}{K}-aware
    model exactly.
  }
  \item{interactions, blocks}{
    Optional \code{\link{interactions}} and \code{\link{blocks}} constraints on
    this forest. The arguments of the same names on \code{\link{dbarts}} are
    the FIRST forest's, so declaring both there is refused as one constraint
    given twice; a second forest's are only expressible here. Holding the two
    forests to different structures is the calibrated-additivity idiom - an
    additive or low-order modulating forest beside a free prognostic one. A
    \code{blocks} partition covers the columns the forest may split on, i.e.
    the \code{vars} subset when one is given.
  }
  \item{amplitude.prior.variance}{
    Prior variance of the \eqn{N(0, \cdot)}{N(0, .)} amplitudes on this
    forest's basis columns, default \code{0.5}. Legal only on a forest given a
    \code{basis}: a forest without one carries a plain scalar amplitude under
    the engine's half-Cauchy scale-mixture prior, whose median is \code{sd}
    rather than a variance. It is a free multiplier on the induced prior: the
    prior standard deviation of the combined location at row \eqn{i}{i} is
    \eqn{\sqrt{\sum_f s_f^2 v_f \|B_f(i,\cdot)\|^2}}{sqrt(sum_f s_f^2 v_f
    ||B_f(i, .)||^2)} over the forests carrying a basis, with \eqn{s_f}{s_f}
    read from \code{$getCalibration(f)[, "prior.scale"]} and \eqn{v_f}{v_f}
    this argument; a basis-free forest's own term is Cauchy and has no
    standard deviation. The budget that sum sits in is set by \code{sd}, whose
    default already divides it among the \eqn{K}{K} forests, so raising this
    argument raises the total rather than redistributing it. Under \code{"probit"} and \code{"logistic"} that
    location IS the latent index and \eqn{\sigma}{sigma} is pinned, so nothing
    in the sampler absorbs a mis-scaled basis; under a gaussian response it is
    in \eqn{\mathrm{sd}(y)}{sd(y)} units and a drawn \eqn{\sigma}{sigma}
    partly does. Every input to that expression is readable off the fitted
    sampler: \eqn{v_f}{v_f} is the \code{amplitude.prior.variance} column of
    \code{$getCalibration(f)} and \eqn{B_f}{B_f} is \code{data@bases[[f]]},
    so the induced prior can be checked against what is in force rather than
    against what the call asked for. See the example below.
  }
  \item{update.amplitude}{
    Whether this forest's amplitudes are redrawn each sweep, default
    \code{TRUE}. \code{FALSE} fixes them at their prior center for the
    sampler's life; the choice is made at creation and cannot be toggled
    afterwards. Declaring it needs a model with amplitudes, so it is refused on
    a single-forest \code{forests}.
  }
}
\details{
  With two forests, the second carrying a two-level factor basis, this is the
  Bayesian causal forest
  \eqn{y = a\,\mu(x) + b_z\,\tau(x) + \epsilon}{y = a mu(x) + b_z tau(x) +
  eps}: a prognostic forest \eqn{\mu}{mu} over every predictor, a modulating
  forest \eqn{\tau}{tau} over the columns \code{vars} allows, and the
  amplitudes \eqn{(a, b_0, b_1)}{(a, b0, b1)} joining them, read back with
  \code{$getForestAmplitudes}. Gaussian, \code{"probit"} and \code{"logistic"}
  responses; under a latent family the combination is the index rather than
  the mean, on the link's own fixed scale, and \code{"aft"}, \code{"ordinal"}
  and \code{"nbinom"} are refused at creation naming what each is missing.

  What the defaults put on the combined location, since under a latent family
  that location IS the index and no drawn \eqn{\sigma}{sigma} stands between it
  and the fitted probabilities. At two forests of the shipped shape - one
  carrying no basis - a probit model's prior puts
  \eqn{P(p < 0.01 \mathrm{~or~} p > 0.99)}{P(p < 0.01 or p > 0.99)} at 0.238,
  which is the shipped single-forest binary default's own 0.239; before the
  \code{sd} defaults above it was 0.376. The \eqn{\sqrt{2/K}}{sqrt(2/K)} factor
  is what holds that as \eqn{K}{K} grows, and it holds it in two different
  senses. When EVERY forest carries a basis the induced prior standard
  deviation of the index is 1.484 latent units at every \eqn{K}{K} - 0.989 of
  the classic \eqn{k = 2}{k = 2} binary node-scale budget - because the whole
  location is then a sum of fixed-variance channels. When one forest carries
  none, its amplitude is Cauchy and has no variance to enter that budget with,
  so the fixed-variance part is BOUNDED by 1.484 rather than pinned at it,
  rising from 0.699 of the budget at \eqn{K = 2}{K = 2} toward 0.989 and never
  reaching it; without the factor it would instead grow past twice the budget
  by ten forests. Read the values in force off
  \code{$getCalibration(f)}'s \code{node.scale.factor} and
  \code{amplitude.prior.scale} columns.

  Both forests' leaf scales come from the model's own calibration map rather
  than from the node prior, which is why a \code{k} hyperprior, a non-default
  \code{k}, and a linear or Gaussian-process node prior are refused when a
  second forest is declared. Every value here is validated at fit time, and
  anything today's engine cannot honour is refused there by name rather than
  dropped.
}
\value{
  A \code{dbartsForest} specification object, resolved when a sampler is built.
}
\seealso{
  \code{\link{dbarts}}, \code{\link{dbartsSpec}}, \code{\link{interactions}},
  \code{\link{blocks}}
}
\examples{
set.seed(0)
n <- 100L
x <- matrix(runif(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
z <- rbinom(n, 1L, 0.5)
y <- 2 * x[, 1] + z * (1 + 2 * x[, 3]) + rnorm(n, 0, 0.2)

sampler <- dbarts(x, y,
                  forests = list(forest(),
                                 forest(basis = ~ factor(z),
                                        vars = c("x1", "x3"),
                                        n.trees = 25L, sd = 1.5)),
                  control = dbartsControl(n.chains = 1L, n.trees = 25L,
                                          n.samples = 20L, n.burn = 20L))
samples <- sampler$run(20L, 20L)
amplitudes <- sampler$getForestAmplitudes()

# the induced prior sd of the combined location, read off the sampler rather
# than recomputed from the call: forest 2 is the one carrying a basis
calibration <- sampler$getCalibration(2L)
basis <- sampler$data@bases[[2L]]
indexSd <- sqrt(calibration[1L, "prior.scale"]^2 *
                calibration[1L, "amplitude.prior.variance"] *
                rowSums(basis^2))
}
\keyword{models}
