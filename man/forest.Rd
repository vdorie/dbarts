\name{forest}
\alias{forest}
\title{Forest Specification for Multi-Forest Models}
\description{
  Build the specification of one forest of a multi-forest model. Pass a list of
  these as the \code{forests} argument of \code{\link{dbarts}} or
  \code{\link{dbartsSpec}}, or write each as a term of a \code{\link{bart}}
  or \code{\link{dbarts}} formula (\sQuote{Details}); either fits the mean as
  a weighted sum of the declared ensembles. Every knob is per forest and carried on the one
  constructor, so the fitting functions grow exactly one argument however many
  forests a model has.

  The constructor is not exported: it resolves by bare name inside the
  arguments that take it, and elsewhere is written
  \code{dbartsForests$forest(...)}; see \code{\link{dbartsForests}}.
}
\usage{
forest(
    vars = NULL, basis = NULL, sd = NULL,
    n.trees = NULL, base = NULL, power = NULL, amplitude = NULL,
    interactions = NULL, blocks = NULL,
    amplitude.prior.variance = NULL)
}
\arguments{
  \item{vars}{
    The predictors this forest splits on, and the one argument that may be
    given without its name: every other argument is given by name. \code{NULL},
    or no first argument at all, leaves the forest every predictor of the fit.

    In a formula given to \code{\link{bart}} or \code{\link{dbarts}} it is
    written as the right-hand side of a model formula of the forest's own:
    \code{forest(x1 + x2)}, \code{forest(log(x1) + factor(g))},
    \code{forest(. - z)}, where \code{.} is every column of \code{data} but the
    response and \code{-} removes a term. The fit's predictors are the terms
    of all of its forests, each once, and a forest splits on the columns of its
    own terms; a term that gives several columns, such as \code{poly(x1, 2)} or
    a factor under \code{factors = "indicators"}, is all of them. Predictors
    may also be named, as \code{forest(c("x1", "x2"))}: a name of a column of
    \code{data} is the term that name written out is, and any other name is a
    column of the predictors as the fit builds them, one of a factor's
    indicator columns for instance. A position is not taken in a formula. An
    \code{offset()} and an intercept term (\code{1}, \code{0} or \code{- 1})
    belong to the fit and are written beside the forests, not inside one. A
    removal written beside the forests takes its term from the forest with no
    \code{basis} and from no other: that forest has the terms of the formula
    with its \code{forest()} replaced by what is inside it and every forest
    with a \code{basis} left out. So
    \code{y ~ forest(.) - x3 + forest(x3, basis = ~ z)} and
    \code{y ~ . - x3 + forest(x3, basis = ~ z)} are one model, in which the
    second forest alone splits on \code{x3}; a removal that leaves the first
    forest no predictor is refused.

    In a \code{forests} list it selects among the fit's predictors, in one of
    two ways. When a name in it is a predictor of the fit, it is read as terms
    in the same way: \code{forest(x1 + x3)}, \code{forest(. - x2)},
    \code{forest(. - log(x1))}. There \code{.} is every column of the fit's
    predictors, and a term names a predictor as the fit holds it, by its term
    label when the fit has a formula and by its column name when it has a
    matrix, written as code or as a backticked name; a term that names no
    predictor is refused, a removed one too. Anything else is a value: it is
    evaluated once, where \code{forest()} is called and at that moment, and
    gives the predictors' names or their positions among the fit's columns:
    \code{forest(c("x1", "x3"))}, \code{forest(c(1, 3))}, or a variable
    holding either. A warning that evaluation raises is kept and raised when a
    fit uses the value. A name or a position given twice is refused, and so
    are a logical, a factor and a formula: the multiplier of a forest is its
    \code{basis}, by name. A column with no name is selected by position or
    by \code{.} only, and a name that two columns share is refused.

    In short, a name of a predictor is that predictor, and anything else is
    read where \code{forest()} is called, at that moment. A forest built in a
    loop, by \code{lapply(names, forest, basis = ~ z)} or by \code{Map()}
    therefore keeps the selection each call was given, whatever its variables
    hold by the time the model is fitted, and a forest saved to a file carries
    its selection with it. A variable that does not exist yet when
    \code{forest()} is called is not looked for later. One case needs care: a
    predictor hides a variable of yours with the same name, so in
    \code{forest(x1)} with a predictor \code{x1} and a variable \code{x1}
    that holds names, the predictor is meant. Where the two could meet, hand
    the value over: \code{do.call(forest, list(vars = nms))}.

    Any forest may be restricted, a single declared forest included: given the
    same residual scale estimate, that fit is the fit on the named columns
    alone, and its other columns report no splits. The estimate is taken from
    every column unless \code{sigest} states it, so the two fits agree draw
    for draw when it does. On a single forest a \code{dart} tree prior beside \code{vars}
    lays its Dirichlet over the named columns, every other column reporting
    probability 0, and \code{split.probs} keeps its ratios among the named
    columns and must give one of them a positive probability (see
    \code{\link{dbartsPriors}}). On the single forest of a hazard fit the
    \code{period} column the fit appends is always allowed, named or not:
    \code{vars} restricts the columns the caller supplied. A hazard fit with
    several forests takes each forest's \code{vars} as written, so a
    restricted forest there splits on \code{period} only when \code{vars}
    names it. The restriction holds for the sampler's
    life; \code{$setModel} refuses a model that states another, or none.
  }
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
    is meaningful: the second level is the one \eqn{b_1}{b1} scales. At
    creation, a column no observation enters - all zeros, or a factor level no
    row takes - is refused, since its amplitude would move under its prior
    alone, as is a basis so large or small that a row's norm is not
    representable; \code{$setForestBasis} keeps an empty factor level (see
    \code{\linkS4class{dbartsSampler}}). A forest
    that declares none takes the implicit intercept its single amplitude
    \eqn{a}{a} scales; every forest past the first needs one, since the
    amplitudes multiplying it are what distinguish it from the first.
    Reaching \code{\link{dbarts}} through an already-built \code{dbartsData}
    \code{formula} is the one route that refuses this argument outright,
    since the declaration would otherwise have nowhere to ride and be
    silently discarded; \code{\link{dbartsSpec}} always takes a pre-built
    data object and installs the declaration instead, replacing whatever
    bases it carried. On the formula interface, an already-evaluated vector or
    matrix must have one row per observation of the FULL \code{data} - the
    same rule a formula's own evaluation already follows - and is then
    restricted to the rows \code{subset} keeps, exactly as a predictor column
    is; a value already restricted to those rows instead, matching the
    subset's count rather than the full data's, is refused by name rather
    than aligned to them by position. With no \code{subset}, or one that
    keeps every row, the two counts agree and nothing changes.

    The right-hand side of a formula is evaluated as R code, and the forest
    does not centre or scale a basis, so a multiplier that should be
    standardized is written so: \code{basis = ~ scale(w)}, or
    \code{basis = ~ cbind(scale(w), scale(v))} for two columns.

    For a \code{forest()} term of a model formula, \code{predict} rebuilds the
    basis at the new rows from what \code{scale()}, \code{poly()},
    \code{ns()} or \code{bs()} computed on the training rows, as
    \code{\link{lm}} does, so a new row is standardized by the training centre
    and scale and not by the new rows' own. This holds for such a call written
    alone, or as an operand of \code{+}, \code{-}, \code{*}, \code{/},
    \code{^}, parentheses or \code{cbind()}; any other expression, such as
    \code{I(w - mean(w))}, \code{scale(w)[, 1]} or \code{abs(scale(w))}, is
    evaluated on the rows given to \code{predict}, as \code{lm} evaluates it,
    so its value at a row depends on the other rows predicted with it; see
    \code{\link[stats]{SafePrediction}}. With \code{subset}, the centre and
    scale are those of the rows \code{subset} kept, where \code{lm} uses all
    of them; rows dropped for a missing response still enter them. A factor
    level the fit never saw is refused at \code{predict}.

    A basis declared through \code{forests =} is not rebuilt: the caller gives
    it at the new rows, and with \code{subset} its formula was evaluated on
    every row of \code{data}.
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
    more specific of the two declarations. \code{\link{bart}} refuses its own
    \code{n.trees} beside one on a formula's forest with no \code{basis}, the
    two being the same count. Where every forest has a \code{basis}, the
    forest written first takes the fitting function's tree count and tree
    prior, so state \code{n.trees} on each forest of such a model.
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
    leaf prior. The two channels default differently, and neither default is a
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
    model exactly. It is one unnamed number, positive and finite, the same for
    every column of the forest's \code{basis}; a vector of any other length,
    a named number and anything that is not a number are refused.
    A model of one forest states its size as the fitting function's
    \code{leaf.prior = normal(sd = )}.

    It is the one argument a live sampler restates:
    \code{$setLeafPrior(forests = list(forest(sd = ), ...))} writes it in the
    channel creation gave the forest, and \code{$getLeafPrior(f)$leaf.prior}
    reads it back as \code{forest(sd = )} (see \code{\link{dbartsSampler}}).
    The channel is fixed at creation: a forest created without a basis keeps
    its half-Cauchy median after \code{$setForestBasis} gives it one, so its
    \code{forest(sd = )} round-trips through \code{$setLeafPrior}, while a
    fresh \code{dbarts()} given the same bases would read the same \code{sd}
    as a leaf-scale factor and build a different prior.
  }
  \item{interactions, blocks}{
    Optional \code{\link{interactions}} and \code{\link{blocks}} constraints on
    this forest. The arguments of the same names on \code{\link{dbarts}} are
    the FIRST forest's, so declaring both there is refused as one constraint
    given twice; a second forest's are only expressible here. Holding the two
    forests to different structures is the calibrated-additivity idiom - an
    additive or low-order modulating forest beside a free prognostic one. A
    \code{blocks} partition covers the columns the forest may split on, i.e.
    the \code{vars} subset when one is given, on the first forest as on any
    other.
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
    read from \code{$getLeafPrior(f)$k.scale} and \eqn{v_f}{v_f}
    this argument; a basis-free forest's own term is Cauchy and has no
    standard deviation. The budget that sum sits in is set by \code{sd}, whose
    default already divides it among the \eqn{K}{K} forests, so raising this
    argument raises the total rather than redistributing it. Under \code{"probit"} and \code{"logistic"} that
    location IS the latent index and \eqn{\sigma}{sigma} is pinned, so nothing
    in the sampler absorbs a mis-scaled basis; under a gaussian response it is
    in \eqn{\mathrm{sd}(y)}{sd(y)} units and a drawn \eqn{\sigma}{sigma}
    partly does. Every input to that expression is readable off the fitted
    sampler: \eqn{v_f}{v_f} is the \code{amplitude.prior.variance} entry of
    \code{$getLeafPrior(f)} and \eqn{B_f}{B_f} is \code{data@bases[[f]]},
    so the induced prior can be checked against what is in force rather than
    against what the call asked for. See the example below.
    A lone forest carrying a \code{basis} is refused, the amplitudes being what
    distinguish a forest from another. For varying coefficients declare an
    intercept forest plus one basis forest per covariate,
    \code{forests = list(forest(), forest(basis = ~ z1), forest(basis = ~ z2))},
    or use one forest with \code{linear()} leaves.
  }
  \item{amplitude}{
    The law of this forest's amplitudes: the one a forest with no
    \code{basis} carries, or one for each column of its \code{basis}. Left out,
    they are drawn each sweep. \code{fixed()} holds them for the sampler's life,
    the choice being made at creation: a forest with no \code{basis} at 1, and
    a forest on a factor, a character or a logical vector at 0 for the first
    level and 1 for each of the others. The forest is then held out of the
    first level and enters at full size for every other level, so it is each
    other level's difference from the first. A basis of several numeric columns
    is held the same way, the first column at 0, so that column is not used. A
    basis of one numeric column would hold the forest at zero and is refused,
    at creation and by \code{$setForestBasis}. \code{fixed} alone is taken as
    \code{fixed()}, and no other value is. Outside the argument that takes it
    the constructor is \code{dbartsPriors$fixed()}. It is stated for a forest
    of a model of several; a model of one forest has no amplitude to hold.
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
  the classic \eqn{k = 2}{k = 2} binary leaf-scale budget - because the whole
  location is then a sum of fixed-variance channels. When one forest carries
  none, its amplitude is Cauchy and has no variance to enter that budget with,
  so the fixed-variance part is BOUNDED by 1.484 rather than pinned at it,
  rising from 0.699 of the budget at \eqn{K = 2}{K = 2} toward 0.989 and never
  reaching it; without the factor it would instead grow past twice the budget
  by ten forests. Read the values in force off
  \code{$getLeafPrior(f)}'s \code{leaf.scale.factor} and
  \code{amplitude.prior.scale} entries.

  Both forests' leaf scales come from the model's own calibration map rather
  than from the leaf prior, which is why a \code{k} hyperprior, a non-default
  \code{k}, a named leaf-prior \code{sd}, and a linear or Gaussian-process
  leaf prior are refused when a second forest is declared. This \code{sd} is
  not the leaf prior's: it states this forest's share of the combined
  location's prior, per unit of basis row norm, where
  \code{normal(sd = )} states a single forest's whole spread. Every value here is validated at fit time, and
  anything today's engine cannot honour is refused there by name rather than
  dropped.

  In a \code{\link{bart}} or \code{\link{dbarts}} formula each
  \code{forest()} term is one forest of the model, without a separate
  \code{forests =} list:
  \code{y ~ forest(x1 + x2) + forest(x1 + x2, basis = ~ z)}. A term stands at
  the top of the formula's right-hand side, joined to the others by \code{+};
  a \code{forest()} crossed with another term, as \code{z:forest(x1 + x2)} or
  \code{z * forest(x1 + x2)}, is refused with the \code{forest()} to write in
  its place. The forest with no \code{basis} is the forest with no multiplier
  and the model's first, wherever it is written. It may be left as plain
  terms, \code{y ~ x1 + x2 + forest(x1 + x2, basis = ~ z)}, which is the same
  model, but a formula has one such forest: plain predictor terms beside a
  \code{forest()} with no \code{basis}, and two such terms, are refused.
  Alone, \code{y ~ forest(x1 + x2)} is the single-forest fit
  \code{y ~ x1 + x2}, in every family that takes a formula. The other forests
  keep the order they are written in. \code{forest()} may also be written
  there as \code{dbartsForests$forest()}. Every argument but the first is
  evaluated in the formula's environment. A
  term's basis is evaluated against the fit's model frame, after
  \code{subset} and any row selection; a \code{basis} declared directly on
  \code{forests =} is instead evaluated against the raw \code{data} and then
  restricted to the same rows - the two routes reach the same rows by
  different timing, so a \code{subset} lines them up either way.
}
\value{
  A \code{dbartsForest} specification object, resolved when a sampler is built.
}
\seealso{
  \code{\link{dbartsForests}},
  \code{\link{dbarts}}, \code{\link{bart}}, \code{\link{dbartsSpec}},
  \code{\link{interactions}}, \code{\link{blocks}}

  \code{\link{bart}}'s \sQuote{Formula Terms} section is the full term
  grammar - the forests of a formula, a forest's predictors, and the refusal
  list - of which the \sQuote{Details} paragraph above is the
  constructor-side summary.
}
\examples{
set.seed(0)
n <- 100L
x <- matrix(runif(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
z <- rbinom(n, 1L, 0.5)
y <- 2 * x[, 1] + z * (1 + 2 * x[, 3]) + rnorm(n, 0, 0.2)

sampler <- dbarts(x, y,
                  forests = list(forest(),
                                 forest(x1 + x3, basis = ~ factor(z),
                                        n.trees = 25L, sd = 1.5)),
                  control = dbartsControl(n.chains = 1L, n.trees = 25L,
                                          n.samples = 20L, n.burn = 20L))
samples <- sampler$run(20L, 20L)
amplitudes <- sampler$getForestAmplitudes()

# the induced prior sd of the combined location, read off the sampler rather
# than recomputed from the call: forest 2 is the one carrying a basis
calibration <- sampler$getLeafPrior(2L)
basis <- sampler$data@bases[[2L]]
indexSd <- sqrt(calibration$k.scale^2 *
                calibration$amplitude.prior.variance *
                rowSums(basis^2))

# a standardized multiplier: predict centres and scales the new rows by the
# training rows' mean and sd of w
d <- data.frame(y = y, x1 = x[, 1], x3 = x[, 3], w = 50 + 10 * rnorm(n))
fit <- bart(y ~ forest(x1 + x3) +
              forest(x1 + x3, basis = ~ scale(w), n.trees = 10L),
            d, n.trees = 10L, n.samples = 10L, n.burn = 10L,
            n.chains = 1L, n.threads = 1L, keepTrees = TRUE)
predict(fit, d[1:3, ])
}
\keyword{models}
