\name{forest}
\alias{forest}
\title{Forest Specification for Multi-Forest Models}
\description{
  One forest of a model of several forests. Such a model fits the mean as a
  sum of forests, each multiplied by columns of the data:
  \code{y ~ forest(x1 + x2) + forest(x1 + x2, basis = dose)} is
  \eqn{a_1 f_1(x) + a_2 \cdot \mathrm{dose} \cdot f_2(x)}{a1 f1(x) + a2 * dose * f2(x)},
  where each \eqn{f}{f} is a forest and each \eqn{a}{a} a coefficient the
  sampler draws. A forest is written as a term of a \code{\link{bart}} or
  \code{\link{dbarts}} formula, as in that example, or as an element of the
  \code{forests} list of \code{\link{dbarts}} or \code{\link{dbartsSpec}}:
  \code{forests = list(forest(), forest(basis = dose))}. Either way
  \code{vars} is the predictors the forest splits on and \code{basis} is what
  it is multiplied by. Every setting is per forest and stated on the forest.

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
    The predictors this forest splits on. In short:
    \itemize{
      \item It is the one argument given without its name:
        \code{forest(x1 + x2)}. Every other argument is given by name.
      \item Left out, the forest splits on every predictor of the fit.
      \item In a formula it is terms joined by \code{+}, as on the right of
        any model formula: \code{forest(x1 + x2)},
        \code{forest(log(x1) + factor(g))}, \code{forest(. - z)}.
      \item In a \code{forests} list the same terms select among the fit's
        predictors, \code{forest(x1 + x3)}, \code{forest(. - x2)}; names or
        positions do too, \code{forest(c("x1", "x3"))},
        \code{forest(c(1, 3))}.
      \item A name of a predictor is that predictor. Anything else is read
        where \code{forest()} is called, at that moment.
    }
    The rest of this item is the detail of each.

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
    \code{basis} and from no forest that names predictors of its own: the
    forest with no \code{basis} has the terms of the formula with its
    \code{forest()} replaced by what is inside it and every forest with a
    \code{basis} left out. So
    \code{y ~ forest(.) - x3 + forest(x3, basis = z)} and
    \code{y ~ . - x3 + forest(x3, basis = z)} are one model, in which the
    second forest alone splits on \code{x3}; a removal that leaves the first
    forest no predictor is refused. Written with no first argument, that
    forest is every predictor the other forests name, less what is removed
    beside it. A removal of a term the forest with no \code{basis} does not
    have is passed by, as R passes it by in a plain formula:
    \code{y ~ x1 + x2 + forest(x3, basis = z) - x3} is the model without the
    removal. \code{\link[stats]{update}} cannot take a term out of a
    \code{forest()}: \code{update(f, . ~ . - x3)} returns \code{f} unchanged
    when \code{x3} stands inside one.

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
    by \code{.} only, and a name that two columns share is refused. A number
    is a position only as the whole argument: among terms, as in
    \code{forest(. - 2)}, it is refused, and a column named \code{2} is
    written backticked. With a factor expanded to indicator columns, removing
    one of them, as \code{forest(. - g.u)}, leaves the forest able to separate
    that level through the others; remove the factor, \code{forest(. - g)}, to
    keep it out.

    A name of a predictor is that predictor, and anything else is
    read where \code{forest()} is called, at that moment. A forest built in a
    loop, by \code{lapply(names, forest, basis = z)} or by \code{Map()}
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
    What the forest is multiplied by. In short:
    \itemize{
      \item A name of a column of \code{data} is that column:
        \code{basis = dose}.
      \item \code{+} separates columns, each with a coefficient of its own:
        \code{dose + age}. \code{I(dose + age)} is one column, their sum.
      \item A factor, a character column or a logical column is one column
        for each of its levels.
      \item Anything that is not a column of \code{data} is read where
        \code{forest()} is called, at that moment, with or without a tilde
        written in place: \code{forest(basis = w)}, \code{forest(basis = ~ w)}.
      \item A formula held in a variable, \code{b <- ~ w} and then
        \code{forest(basis = b)}, is read the same way, when \code{forest()}
        is called, and is left as you made it.
      \item \code{cbind()} is refused: write the columns with \code{+}.
    }
    The rest of this item is the detail of each.

    \strong{What it does.} With a basis of columns
    \eqn{B_1, \ldots, B_K}{B_1, ..., B_K} the forest \eqn{f}{f} enters the fit
    as \eqn{(\sum_k a_k B_k(x_i)) f(x_i)}{(sum_k a_k B_k(x_i)) f(x_i)}, each
    column with a coefficient \eqn{a_k}{a_k} of its own, called its
    amplitude. A forest with no basis is multiplied by a single coefficient.
    A model has one forest with no basis, or none; every other forest states
    one, and a model whose only forest has a basis is refused.

    \strong{Writing it.} A basis is the right-hand side of a model formula,
    without the tilde, and means what it means in \code{\link{lm}}:
    \itemize{
      \item a name is a column: \code{basis = dose};
      \item \code{+} separates columns: \code{dose + age} is two columns;
      \item arithmetic is written inside \code{I()}: \code{I(dose + age)} is
        one column, the sum, and \code{I(dose / 30)} is \code{dose}
        rescaled;
      \item a function of columns is a term, as in \code{lm}:
        \code{log(dose)}, \code{scale(age)}, and \code{poly(dose, 2)},
        which is two columns;
      \item a factor, a character column or a logical column gives one column
        for each of its levels, none left out: \code{factor(z)} for a
        \code{z} of 0 and 1 is the two columns whose coefficients are
        \eqn{(b_0, b_1)}{(b0, b1)}, in the order of the levels;
      \item \code{dose:age} is one column, the product;
      \item there is no constant column unless \code{1} asks for one:
        \code{1 + dose} is a column of ones and \code{dose}, while
        \code{0 + dose} and \code{dose - 1} are \code{dose}.
    }
    A basis is either one factor (or character or logical column) or numeric
    columns; the two are not mixed in one forest. The forest does not centre
    or scale a basis, so a multiplier that should be standardized is written
    so: \code{basis = scale(age)}.

    Some forms a model formula would take are refused, each by name and with
    what to write instead: \code{*} between terms (write \code{I(dose * age)}
    for the product, or \code{dose + age + I(dose * age)} for all three
    columns); a term multiplied or divided by a number (write
    \code{I(dose / 30)}); \code{-} between terms, other than \code{- 1};
    \code{^}, \code{/}, \code{\%in\%} and \code{|} between terms; \code{.};
    an \code{offset()}; \code{cbind()} (write \code{+}); and a term that
    calls one of \code{normal}, \code{fixed}, \code{student}, \code{cauchy},
    \code{linear}, \code{gp}, \code{cgm}, \code{dart}, \code{chisq},
    \code{chi}, \code{invchi}, \code{forest} or \code{varianceForest}, since
    a prior is not stated on a term of a basis.

    \strong{Names of the columns.} They are named as
    \code{coef(lm(y ~ 0 + <basis>))} names them: \code{dose}, \code{age};
    \code{poly(dose, 2)1}, \code{poly(dose, 2)2}; \code{factor(z)0},
    \code{factor(z)1}; \code{(Intercept)} for the column of ones. They are the
    column names of \code{data@bases[[f]]} on the sampler. A forest's label,
    where its place in a \code{forests} list has no name, is the text of its
    basis, \code{"dose + age"}, or \code{forest<i>} for a forest with no basis
    or with one given as a value; see \code{\link{dbarts}}.

    \strong{Rows.} A basis has a value for every row of \code{data}, and is
    cut to the rows the fit keeps: \code{data} and \code{subset} are each
    evaluated once, for the fit, and the basis is given the same rows in the
    same order, less those the \code{na.action} drops. So with
    \code{data = d[sample(nrow(d)), ]} or \code{subset = sample(n, 100)} the
    basis has the rows the fit has. \code{scale()}, \code{poly()} and anything
    else a term computes across rows, \code{I(age - mean(age))} among it, are
    computed on every row of \code{data}, as in \code{lm}, whatever
    \code{subset} keeps. A level of a factor that no kept row has gives no
    column. Refused are a missing value in a row the fit keeps, where
    \code{lm} would drop the row; a factor, a character column or a logical
    column with one level on the rows kept; and a numeric column that is zero
    on all of them.

    \strong{New rows.} \code{predict} builds the basis again from the new
    rows, as \code{predict} does for \code{lm}. Taken from the fitted rows
    are the centre, scale and knots of \code{scale()}, \code{poly()},
    \code{ns()} and \code{bs()}, and the levels of a factor. Everything else
    is computed again on the rows given: \code{I(age - mean(age))} uses the
    mean of the new rows, and \code{factor(dose > median(dose))} their
    median; see
    \code{\link[stats]{SafePrediction}}. Every column of \code{data} the basis
    names must be among the new rows, a factor's columns keep the order of the
    fit's levels, and a level the fit did not have is refused. A basis given
    as a value cannot be built again: give it at the new rows through
    \code{predict}'s \code{bases} argument.

    \strong{Where a name is found.} It depends on who wrote the code.

    In a call of \code{forest()} that you write, as in a \code{forests} list
    or for a forest built ahead of the fit, a name of a column of \code{data}
    is that column. Every other variable of yours that the code uses is
    copied when \code{forest()} is called. This is so whether or not a tilde
    is written in place: \code{forest(basis = ~ dose + age)} is
    \code{forest(basis = dose + age)}. A forest built in a loop therefore
    keeps what each call was given,
    \preformatted{
    fs <- list(forest())
    for (k in c(10, 30))
      fs[[length(fs) + 1]] <- forest(x1, basis = I(dose / k))
    }
    is one forest multiplied by \code{dose / 10} and one by \code{dose / 30},
    at the fit and at \code{predict}, whatever \code{k} holds later, and a
    saved forest carries its numbers with it. A variable that does not exist
    when \code{forest()} is called is not looked for later, and is refused.
    What R or an attached package supplies, a function such as \code{scale}
    or a constant such as \code{pi}, is not copied: a basis that names a
    column of \code{data} looks it up when the model is fitted. With no data
    frame, as with the matrix interface and
    \code{\link{dbartsSpec}}, there is no column to name and every variable is
    yours.

    One case needs care: a column of \code{data} hides a variable of yours
    with the same name, whether that variable holds a column, a set of names
    or a one-sided formula. Where the two could meet, hand the value over:
    \code{do.call(forest, list(basis = w))} for a column you hold,
    \code{list(basis = as.name(nm))} for a column of \code{data} by its name,
    \code{list(basis = f)} for a one-sided formula, and
    \code{list(vars = nms)} for predictors by name.

    A formula made elsewhere, held in a variable or handed over as above, is
    read in the same way and at the same moment, when \code{forest()} is
    called: a name of a column of \code{data} is that column, and every other
    variable the formula uses is copied, as it is then, from where the
    formula was made. What counts is when \code{forest()} is called, not when
    the formula was made. In
    \preformatted{
    fs <- list(forest())
    for (k in c(10, 30)) {
      f <- ~ I(dose / k)
      fs[[length(fs) + 1]] <- forest(x1, basis = f)
    }
    }
    each forest has its own \code{k}, because \code{forest()} is called
    inside the loop. Formulas made in a loop and given to \code{forest()}
    only after it all see the last \code{k}: a formula holds its environment
    and no value, as everywhere in R. A variable changed after
    \code{forest()} is called changes neither the fit nor a prediction. The
    formula itself is not touched: it is the formula you made, in its
    environment, and nothing is assigned there.

    A \code{forest()} term inside the formula of \code{\link{bart}} or
    \code{\link{dbarts}} is part of that formula and is read with it, as R
    reads a formula, with a tilde or without: in \code{data} and then where
    the formula was written, when the model is fitted and again by
    \code{predict}. Change a number such a term uses between the fit and
    the prediction and the prediction changes, with no message.

    An argument that \code{lapply()}, \code{Map()} or another function of
    base R writes for you, as in \code{lapply(columns, forest, vars = "x1")},
    is a value handed over, whatever the data's columns are called. A name
    or a call that such a function is handed and passes on, as in
    \code{Map(forest, nms, MoreArgs = list(basis = quote(dose)))}, is code,
    as it is through \code{do.call()}.

    \strong{A value.} What is handed over as an object, through
    \code{do.call()} or in a call built by a program, is used as it is: a
    numeric vector or matrix is its columns, and a factor, a character vector
    or a logical vector is its levels. Its columns keep their names when every
    column has one, and have none when any column lacks one; two columns of
    one name are refused. Its forest's label is \code{forest<i>}. It has one
    row for every row of \code{data} and is cut to the rows the fit keeps; a
    value already cut to those rows is refused by name. A level that no row
    of \code{data} has is refused, while one that only \code{subset} empties
    keeps a column of zeros, a value not being built again on the rows kept.
    A value that the rows kept leave with one level, or one numeric column
    that is zero on all of them, is refused as the same basis written as code
    is.
    A single string or number is refused: \code{basis = "dose"} names no
    column, so write \code{basis = dose}.

    A basis is stated here or on the data object, not both:
    \code{\link{dbarts}} given an already-built \code{dbartsData} refuses
    this argument, the data object's own \code{bases} being where it is
    stated then, while \code{\link{dbartsSpec}} takes a data object and
    installs the bases declared here in place of those it carried. After
    creation a basis is replaced by \code{$setForestBasis}, whose formula is
    read at once, where it was made; see \code{\linkS4class{dbartsSampler}}.
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
    model exactly. It is one unnamed number, positive and finite. For a basis
    of several columns that one number is the size stated for each of them,
    so \code{forest(basis = dose + age, sd = 2)} states 2 for the column
    \code{dose} and 2 for the column \code{age}; to size two columns
    differently, rescale one of them in the basis, as \code{I(dose / 30)}. A
    vector of any other length, a named number and anything that is not a
    number are refused.
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
    \code{forests = list(forest(), forest(basis = z1), forest(basis = z2))},
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
  \code{y ~ forest(x1 + x2) + forest(x1 + x2, basis = z)}. A term stands at
  the top of the formula's right-hand side, joined to the others by \code{+};
  a \code{forest()} crossed with another term, as \code{z:forest(x1 + x2)} or
  \code{z * forest(x1 + x2)}, is refused with the \code{forest()} to write in
  its place. The forest with no \code{basis} is the forest with no multiplier
  and the model's first, wherever it is written. It may be left as plain
  terms, \code{y ~ x1 + x2 + forest(x1 + x2, basis = z)}, which is the same
  model, but a formula has one such forest: plain predictor terms beside a
  \code{forest()} with no \code{basis}, and two such terms, are refused.
  Alone, \code{y ~ forest(x1 + x2)} is the single-forest fit
  \code{y ~ x1 + x2}, in every family that takes a formula. The other forests
  keep the order they are written in. \code{forest()} may also be written
  there as \code{dbartsForests$forest()}. A term's predictors and its basis
  are read as code; every other argument is evaluated where the formula was
  written. The same basis written in a formula's term and in a \code{forests}
  list is the same columns on the same rows, and so the same model.
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
                                 forest(x1 + x3, basis = factor(z),
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

# the same model written in a formula, the first forest written out; the
# columns of a basis are named as lm() would name them
d <- data.frame(y = y, x1 = x[, 1], x2 = x[, 2], x3 = x[, 3], z = z,
                w = 50 + 10 * rnorm(n))
written <- dbarts(y ~ forest(x1 + x2 + x3) +
                    forest(x1 + x3, basis = factor(z), n.trees = 25L, sd = 1.5),
                  d, control = dbartsControl(n.chains = 1L, n.trees = 25L))
colnames(written$data@bases[[2L]])

# a standardized multiplier: predict centres and scales the new rows by the
# training rows' mean and sd of w
fit <- bart(y ~ forest(x1 + x3) +
              forest(x1 + x3, basis = scale(w), n.trees = 10L),
            d, n.trees = 10L, n.samples = 10L, n.burn = 10L,
            n.chains = 1L, n.threads = 1L, keepTrees = TRUE)
predict(fit, d[1:3, ])

# forests built in a loop: each keeps the number its own call was given
forest <- dbartsForests$forest
forests <- list(forest())
for (k in c(25, 100))
  forests[[length(forests) + 1L]] <- forest(x1, basis = I(w / k))
looped <- dbarts(y ~ x1 + x3, d, forests = forests,
                 control = dbartsControl(n.chains = 1L, n.trees = 10L))
vapply(looped$data@bases[-1L], function(basis) basis[1L, 1L], 0) * c(25, 100)

# two columns, a coefficient each; I() holds arithmetic
twoColumns <- dbarts(y ~ forest(x1 + x3) + forest(x1, basis = z + I(w / 50)),
                     d, control = dbartsControl(n.chains = 1L, n.trees = 10L))
colnames(twoColumns$data@bases[[2L]])
}
\keyword{models}
