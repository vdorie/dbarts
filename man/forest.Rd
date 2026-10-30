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
    proposal.probs = NULL, n.perturb.cuts = NULL, treeShift = NULL)
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
    \code{\link{dbartsPriors}}). On every forest of a hazard fit the
    \code{period} column the fit appends is always allowed, named or not:
    \code{vars} restricts the columns the caller supplied. The restriction
    holds for the sampler's
    life; \code{$setModel} refuses a model that states another, or none.
  }
  \item{basis}{
    What the forest is multiplied by. In short:
    \itemize{
      \item A name of a column of \code{data} is that column:
        \code{basis = dose}.
      \item A basis is one numeric column, or a factor, a character column or
        a logical column of two levels, with a coefficient for each level. A
        number is one column, whatever its values: \code{basis = z} for a
        \code{z} of 0s and 1s is one coefficient, and
        \code{basis = factor(z)} is the two levels.
      \item \code{I(dose + age)} is one column, their sum. Several columns,
        \code{dose + age}, and three or more levels are refused: give each
        column, or each level but the first, a forest of its own.
      \item Anything that is not a column of \code{data} is read where
        \code{forest()} is called, at that moment, with or without a tilde
        written in place: \code{forest(basis = w)}, \code{forest(basis = ~ w)}.
      \item A formula held in a variable, \code{b <- ~ w} and then
        \code{forest(basis = b)}, is read the same way, when \code{forest()}
        is called, and is left as you made it.
      \item \code{cbind()} is refused: write its one column alone, or give
        each column a forest of its own.
    }
    The rest of this item is the detail of each.

    \strong{What it does.} With a basis of columns
    \eqn{B_1, \ldots, B_K}{B_1, ..., B_K} the forest \eqn{f}{f} enters the fit
    as \eqn{(\sum_k a_k B_k(x_i)) f(x_i)}{(sum_k a_k B_k(x_i)) f(x_i)}, each
    column with a coefficient \eqn{a_k}{a_k} of its own, called its
    amplitude: \eqn{K = 1}{K = 1} for a numeric column, and the two columns
    of a factor's levels, \eqn{K = 2}{K = 2}. A forest with no basis is
    multiplied by a single coefficient.
    A model has one forest with no basis, or none; every other forest states
    one, and a model whose only forest has a basis is refused.

    \strong{Writing it.} A basis is the right-hand side of a model formula,
    without the tilde, and means what it means in \code{\link{lm}}:
    \itemize{
      \item a name is a column: \code{basis = dose};
      \item \code{+} separates columns: \code{dose + age} is two columns,
        which this version refuses (below);
      \item arithmetic is written inside \code{I()}: \code{I(dose + age)} is
        one column, the sum, and \code{I(dose / 30)} is \code{dose}
        rescaled;
      \item a function of columns is a term, as in \code{lm}:
        \code{log(dose)}, \code{scale(age)} and \code{poly(dose, 1)};
      \item a factor, a character column or a logical column of two levels
        gives one column for each, none left out: \code{factor(z)} for a
        \code{z} of 0 and 1 is the two columns whose coefficients are
        \eqn{(b_0, b_1)}{(b0, b1)}, in the order of the levels;
      \item \code{dose:age} is one column, the product;
      \item there is no constant column unless \code{1} asks for one, and
        \code{1 + dose}, a column of ones beside \code{dose}, is refused,
        while \code{0 + dose} and \code{dose - 1} are \code{dose}.
    }
    A basis is either one factor (or character or logical column) or numeric
    columns; the two are not mixed in one forest. The forest does not center
    or scale a basis, so a multiplier that should be standardized is written
    so: \code{basis = scale(age)}. A numeric column constant over the
    observations in the likelihood has no standard deviation to take the
    default \code{sd} from, and is refused unless \code{sd} is stated.

    \strong{One size.} In this version a basis is one numeric column, or a
    factor, a character or a logical vector of two levels, counted on the
    rows the fit keeps once a level they leave empty is dropped. Three or more
    levels are refused, with the form to write, a forest for each level but
    the first: \code{forest(x1, basis = arm == "b") + forest(x1, basis = arm
    == "c")}. Several numeric columns are refused however they are written,
    \code{dose + age}, \code{poly(dose, 2)}, \code{splines::ns(dose, 3)}, a
    matrix, or a factor's indicator columns handed over as numbers, with the
    form to write, a forest for each column:
    \code{forest(x1, basis = dose) + forest(x1, basis = age)}; indicator
    columns are handed over as the factor itself. \code{1 + dose} is refused
    for its constant column: write \code{basis = dose}. Each is refused drawn
    or held, at every place a basis is given, and a basis handed over as a
    value is told the form in words.

    \strong{Which kind it is.} The class of the basis decides whether the
    forest multiplies levels or numbers, the same way at every place a basis
    is given: \code{forest()} written as code or with a value,
    \code{\link{dbartsData}}'s \code{bases}, \code{$setForestBasis} and
    \code{predict}'s \code{bases}. A factor, ordered or not, a character
    vector and a logical vector are levels, with a coefficient for each: in
    the order of a factor's levels, in sorted order for a character vector,
    \code{FALSE} then \code{TRUE} for a logical vector. An integer or a
    double is numbers, whatever its values, and a numeric matrix is its
    columns. The values are never looked at to decide: a column of 0s and 1s
    is one column with one coefficient, and is two levels only when written
    \code{factor(z)} or as a logical vector. A matrix of characters or of
    logicals, and a factor matrix, is refused, one of a single column
    included, since a vector is levels and a matrix is numbers: hand over its
    column as a vector, \code{cm[, 1]}. An array of one dimension is its
    vector. A data frame and any other
    class are refused too. Handed over as a value, the refusal names the
    class; written as code, a class a model frame does not take, a data frame
    or a list, is refused in R's own words.

    The kind is fixed when the forest is created. A forest created with no
    basis cannot be given one, and a forest of levels keeps the same levels,
    in the same order, for as long as the sampler lives, so that coefficient
    \eqn{j}{j} multiplies level \eqn{j}{j}'s rows after any replacement. See
    \code{$setForestBasis} in \code{\linkS4class{dbartsSampler}} and the
    \code{bases} argument of \code{predict} in \code{\link{bartBT}} for what
    a replacement takes.

    A level that no row the fit keeps has gives no column and no message,
    written as code or handed over as a value: \code{factor(z, levels = 0:2)}
    of a \code{z} with no 2, or a factor of three levels under a
    \code{subset} that leaves two of them, is the two columns of the levels
    that remain. The count is made after the rows are cut, and fewer than two
    levels left is refused. A row of weight zero is a row the fit keeps.

    A missing value is looked for on the rows the fit keeps, after
    \code{subset} and the \code{na.action}, as \code{lm} looks for one in its
    variables. On a kept row it is refused. On a row the fit drops it is not
    looked at, so a basis fits the same way written as code or handed over as
    a value. \code{NA} is not a level: a factor made with \code{addNA()} or
    \code{exclude = NULL} that has a kept row at \code{NA} is refused, and
    the remedy is to give that level a name, as in
    \code{levels(f)[is.na(levels(f))] <- "none"}. With no kept row at it the
    \code{NA} level is dropped like any other level no row has.

    Some forms a model formula would take are refused, each by name and with
    what to write instead: \code{*} between terms (write \code{I(dose * age)}
    for the product, or give \code{dose}, \code{age} and
    \code{I(dose * age)} a forest each for all three); a term multiplied or divided by a number (write
    \code{I(dose / 30)}); \code{-} between terms, other than \code{- 1};
    \code{^}, \code{/}, \code{\%in\%} and \code{|} between terms; \code{.};
    an \code{offset()}; \code{cbind()} (write its one column alone, or give
    each column a forest of its own); and a term that
    calls one of \code{normal}, \code{fixed}, \code{student}, \code{cauchy},
    \code{linear}, \code{gp}, \code{cgm}, \code{dart}, \code{chisq},
    \code{chi}, \code{invchi}, \code{forest} or \code{varianceForest}, since
    a prior is not stated on a term of a basis.

    \strong{Names of the columns.} They are named as
    \code{coef(lm(y ~ 0 + <basis>))} names them: \code{dose};
    \code{factor(z)0}, \code{factor(z)1}; and, for a basis of several
    columns, which is read and named before it is refused and whose refusal
    quotes the names, \code{dose}, \code{age}, \code{poly(dose, 2)1},
    \code{poly(dose, 2)2} and \code{(Intercept)} for the column of ones. They
    are the column names of \code{data@bases[[f]]} on the sampler. A forest's label,
    where its place in a \code{forests} list has no name, is the text of its
    basis, \code{"scale(age)"}, or \code{forest<i>} for a forest with no basis
    or with one given as a value; see \code{\link{dbarts}}. Wherever a method
    or an \code{extract} or \code{predict} takes a \code{forest}, a string is
    a forest's label (failing an exact match, the label that is the same code,
    so \code{"I(dose / 30)"} finds \code{I(dose/30)}), a number is its
    position as base R reads a list element by number, and \code{forest<i>} names position
    \emph{i}, though a per-forest margin of a fit is named by the forests'
    labels. A string is never a
    position: \code{"2"} is a forest labeled \code{"2"} or none. A string
    that is one forest's label and another position's name, or the same code
    as two labels, is refused.

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
    are the center, scale and knots of \code{scale()}, \code{poly()},
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
    is written in place: \code{forest(basis = ~ I(dose + age))} is
    \code{forest(basis = I(dose + age))}. A forest built in a loop therefore
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
    or a logical vector is its levels, as written as code. Its column keeps
    the name it has. Its forest's label is \code{forest<i>}. It has one
    row for every row of \code{data} and is cut to the rows the fit keeps; a
    value already cut to those rows is refused by name. A level that no row
    the fit keeps has gives no column, whether \code{data}, \code{subset} or
    the \code{na.action} is what empties it. A value that the rows kept leave with one level, or one numeric column
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

    A lone forest carrying a \code{basis} is refused, the coefficients being
    what distinguish a forest from another. For varying coefficients declare
    an intercept forest plus one basis forest per covariate,
    \code{forests = list(forest(), forest(basis = z1), forest(basis = z2))},
    or use one forest with \code{linear()} leaves.
  }
  \item{n.trees, base, power}{
    This forest's tree count and tree-structure prior. Each one stated here
    governs this forest, and each one left out takes the default of the
    forest's kind, whatever the forest's place. A forest with a \code{basis}
    has 50 trees at \code{base = 0.25}, \code{power = 3} - shallower and
    fewer than a prognostic forest's, the modulating surface normally being
    the smoother of the two. The forest with no \code{basis} has the fitting
    function's: \code{\link{bart}}'s \code{n.trees} or the control's, and
    \code{tree.prior}. A value stated here on that forest governs it over the
    control's \code{n.trees} and the tree prior's \code{base} or
    \code{power}, being the more specific of the two declarations, while
    \code{\link{bart}} refuses its own \code{n.trees} beside one stated here,
    the two being the same count. Where every forest has a \code{basis} the
    fitting function's tree count and tree prior have no forest to belong to,
    and one that is stated is refused, the default's own value included:
    state it on a forest. A sampler's \code{control@n.trees} and
    \code{model@tree.prior} hold its first forest's.
  }
  \item{sd}{
    The size of what this forest contributes, in the response's units: the
    units of \eqn{y}{y} for a continuous response, and of the latent index
    under \code{"probit"} and \code{"logistic"}, where the link's own error
    has standard deviation 1 and \eqn{\pi/\sqrt{3}}{pi/sqrt(3)}. One number;
    a name on it is ignored. What it is the size of depends on what
    multiplies the forest.

    \strong{No basis.} The forest enters as \eqn{a\,f(x)}{a f(x)}, and
    \code{sd} is the size of \eqn{a\,f(x)}{a f(x)}.

    \strong{A factor} (a factor, a character vector or a logical vector, of
    two levels). The forest enters with one coefficient for each level, and
    \code{sd} is the size of the difference between the two: the treatment
    effect where the levels are control and treated.

    \strong{A number} \eqn{w}{w}. A column that is an ordinary covariate
    should be centered: \code{basis = scale(w)}, or
    \code{basis = I(w - c)} for a reference value \eqn{c}{c}. The other
    forests describe the response where \eqn{w}{w} is 0, and a column that
    never comes near 0, an age or a year, makes them an extrapolation and the
    fit poor, with nothing in the chain to show it. A dose or an exposure,
    where 0 means none, is left as it is. The forest enters as
    \eqn{w\,(a\,f(x))}{w (a f(x))}, and \code{sd} is the size of
    \eqn{a\,f(x)}{a f(x)}, the change in the response for one unit of
    \eqn{w}{w}; the term at an observation has size \eqn{|w|}{|w|}
    \code{sd}. The column is used as given: it is never centered and never
    divided, so a stated \code{sd} is per unit of \eqn{w}{w}, and the same
    \code{sd} with \eqn{w}{w} in grams and in milligrams states two priors a
    thousand times apart.

    Where the coefficient is held (\code{amplitude = fixed()}) it is 1, what
    it multiplies is the basis, and "size" is the prior standard deviation,
    exactly. Where it is drawn the size is itself uncertain, and \code{sd} is
    its prior median.

    Not stated, \code{sd} for a forest with no basis is twice
    \eqn{\mathrm{sd}(y)}{sd(y)} under a continuous response where its
    coefficient is drawn, the size \code{\link{bart}} gives its forest (a
    quarter of the response's range) where it is held, and one unit of the
    link's error under \code{"probit"} and \code{"logistic"} either way;
    \eqn{\mathrm{sd}(y)}{sd(y)}, or that unit, for a factor; and for a number
    that same size for one standard deviation of the column,
    \eqn{\mathrm{sd}(y) / \mathrm{sd}(w)}{sd(y) / sd(w)} per unit, held or
    drawn. None of these depends on how many forests the model has. A default
    fit is then the same fit whatever units \eqn{w}{w} is in.
    \eqn{\mathrm{sd}(y)}{sd(y)} and \eqn{\mathrm{sd}(w)}{sd(w)} are taken
    once, over the observations in the likelihood when the sampler is created
    (those \code{subset} and the missing-value handling keep, with a positive
    weight, and not masked out), and kept: replacing the column, the weights
    or the mask later does not move the prior, unless the column is replaced
    with \code{$setForestBasis(updateBasisScale = TRUE)}, and a copy, a saved
    and reloaded sampler and a sampler made again from its own control, model
    and data keep it. Where those observations hold fewer than two distinct
    values of the response, both are taken over every observation, under one
    warning (two by the route of \code{\link{dbartsSpec}} and \code{new()},
    see its \code{active}). \code{"probit"} and \code{"logistic"} read no
    \eqn{\mathrm{sd}(y)}{sd(y)}, and where fewer than two observations are
    in, \eqn{\mathrm{sd}(w)}{sd(w)} is taken over every one, under one
    warning; two or more that hold one value of \eqn{w}{w} are refused. \eqn{\mathrm{sd}(y)}{sd(y)} is of the response net of its
    offset; on a gaussian response that does not vary over those observations
    it is a quarter of the window of width one the response is then scaled
    by. \eqn{\mathrm{sd}(w)}{sd(w)} is kept with the data, in its
    \code{basis.scale} slot (see \code{\link{dbartsData}}), and
    \code{$getLeafPrior(f)$basis.scale} reports it. A caller who replaces a
    column during burn-in and wants the default of the new one replaces it
    with \code{$setForestBasis(f, w.new, updateBasisScale = TRUE)}, which takes
    \eqn{\mathrm{sd}(w)}{sd(w)} again over the observations in the
    likelihood then (see \code{\linkS4class{dbartsSampler}}). A column that
    is constant over those observations has no standard deviation and needs
    \code{sd} stated.

    A treatment belongs in a factor or a logical vector, where \code{sd} is
    the treatment effect. Written as a 0/1 number it is scaled as any number:
    \eqn{\mathrm{sd}(z)}{sd(z)} is 0.5, 0.4 and 0.22 at 50, 20 and 5 percent
    treated, so the default effect is 2, 2.5 and 4.6 times the one
    \code{factor(z)} gets; \code{sd = sd(y)} on the number states the
    factor's. The same holds for a coding written as a number and held,
    \code{basis = z - 0.5} with \code{amplitude = fixed()}.

    A held coefficient cannot adapt to an \code{sd} that is too small. At the
    default it is sized for a slope of about one \eqn{\mathrm{sd}(y)}{sd(y)}
    per \eqn{\mathrm{sd}(w)}{sd(w)}; where the effect may be larger, under
    \code{"probit"} and \code{"logistic"} above all, state \code{sd} or let
    the coefficient be drawn.

    A default is not a stated number: \code{$getLeafPrior(f)$sd} reports the
    number of the response's units it comes to, a numeric forest's named by
    its column, and \code{$getLeafPrior(f)$sd.stated} is \code{FALSE}. A
    stated number is kept as written and divided by the unit once, at each
    creation. A control or a \code{forests} list taken to other data states
    the same numbers in that data's units, and the sampler made from it
    derives its own. It is one number, positive and finite. A vector of any
    other length and anything that is not a number are refused. The one
    forest of a model of one forest takes it as the fitting function's
    \code{leaf.prior = normal(sd = )}, in the same units, and stating the size
    at both places is refused.

    It is the one argument a live sampler restates:
    \code{$setLeafPrior(forests = list(forest(sd = ), ...))} states it from
    then on, also where it is the number in force, and leaves the recorded
    \eqn{\mathrm{sd}(w)}{sd(w)}. \code{$getLeafPrior(f)$leaf.prior} returns
    the statement, \code{forest(sd = )} where one was stated and
    \code{forest()} where none was, and \code{forest()} written back leaves
    its forest as it is (see \code{\link{dbartsSampler}}). A forest's kind is
    fixed at creation: \code{$setForestBasis} refuses to give a basis to a
    forest created without one, and moves no forest's prior.

    \code{print} of a fit and \code{show} of a sampler give a short block of
    one line for each forest of a model of several: its label, what it
    multiplies, for a numeric column the mean and sd it has over the
    observations in the likelihood as they stand, and the sd the forest
    resolved to. For a numeric column that is its sd per sd of the column as
    it now stands, named by the forest's label, the size of the forest's term
    over one sd of the column in the response's units: for a default it is
    \eqn{\mathrm{sd}(y)}{sd(y)} (the link's error sd under \code{"probit"} and
    \code{"logistic"}) wherever the column, and the observations it is read
    over, are those the default was taken from, and for a stated \code{sd} it
    is that \code{sd} times the column's sd. Where the column is constant it is the \code{sd} per
    unit. For a factor it is the sd between the two levels. The sd per unit
    of a numeric column, which is what \code{forest(sd = )} states, is in the
    full block and in the reader. The last line of the short block says where
    each forest's prior is to be had in full: \code{summary()} of a fit,
    which prints the full block, or \code{$getLeafPrior()} of a sampler,
    which returns it. The
    \code{summary} of a fit and the creation of a sampler with
    \code{verbose = TRUE} give the full block. Its header gives the unit the
    sds are in: \eqn{\mathrm{sd}(y)}{sd(y)}, as taken above, or the link's
    error sd. A numeric line gives the sd per unit and per sd of the column
    and adds the column's range; every line adds whether the sd was stated or
    is the default, and whether the coefficient is held. A default whose
    column has another sd over the observations now in the likelihood than
    the one it was taken from (the column replaced, the weights or the active
    rows changed, or, under \code{"gaussian"}, every observation read because
    those in the likelihood held fewer than two values of the response) says
    so. A fit made with \code{keepFits = FALSE} keeps no bases or labels and
    gives no block.
    Two things to look for. The mean of a numeric column, in either block, and
    its range, in the full one: a mean far from 0 is a column that never comes
    near 0, an age or a year, and should be centered, unless 0 means none, as
    for a dose; and one wild value shows in the range, and moves the mean and
    sd. And the figure per sd of the column, in either block, the size of the
    forest's term over one sd of the column (\code{sd} times the column's
    sd): for a default it is \eqn{\mathrm{sd}(y)}{sd(y)} on every numeric line
    whose column, and the observations it is read over, are those the default
    was taken from, so a line that differs has had its column or its
    observations changed, and for a stated \code{sd} it is where a number in
    the wrong units shows. Nothing in the fit warns of either.
  }
  \item{interactions, blocks}{
    Optional \code{\link{interactions}} and \code{\link{blocks}} constraints on
    this forest. The arguments of the same names on \code{\link{dbarts}} are
    those of the forest with no \code{basis}, wherever it stands, so one
    stated there and on that forest is refused as one constraint given
    twice. A forest with a \code{basis} states its own here, and where every
    forest has one those arguments of the fitting function are refused.
    Holding the two forests to different structures is the
    calibrated-additivity idiom - an additive or low-order modulating forest
    beside a free prognostic one. A \code{blocks} partition covers the columns
    the forest may split on, i.e. the \code{vars} subset when one is given,
    and the forest's own tree count.
  }
  \item{amplitude}{
    The prior on this forest's coefficient: the one a forest with no
    \code{basis} carries, or one for each level of its \code{basis}. Left out,
    it is drawn each sweep, under a normal prior on a forest with a
    \code{basis} and, on a forest with none, a half-Cauchy prior whose median
    is the forest's \code{sd}. \code{student(df)} draws it under a Student-t
    prior with \code{df} degrees of freedom, finite and at least 1;
    \code{student()} is \code{student(1)}, the half-Cauchy. With a basis the
    t is scaled so that its median size is the normal's, so the forest's
    \code{sd} keeps its meaning, the prior median of what the forest
    contributes, and only the tails move: a coefficient ten times its median
    size has prior probability 0.06 under \code{student(1)}, 0.005 under
    \code{student(3)} and about \eqn{10^{-11}}{1e-11} under the normal. That
    tail lets the coefficient follow a multiplier whose scale is far from what
    the \code{sd} states, at the cost of slower mixing along the ridge the
    coefficient shares with the forest's leaves. With no basis the
    coefficient's size is the forest's size, its median the \code{sd} or, where
    none is stated, the size the forest has by default, so
    \code{student(1)} restates that default and Bayesian causal forests'
    prognostic prior is written \code{amplitude = student(1), sd = 2 * sd(y)}.
    The draw is exact, a normal coefficient given a variance drawn from its
    inverse-gamma conditional. \code{fixed()} holds it for the sampler's life,
    the choice being made at creation, wherever the forest stands: a held
    coefficient is 1, and what it multiplies is the basis, so a forest with no
    \code{basis} enters at 1, a forest on a numeric column at the column as
    given, and a forest on a factor, a character or a logical vector of two
    levels at 0 for the first level and 1 for the second. The forest is then held out of the first level and enters at full
    size for the second, so it is the second level's difference from the
    first. Another coding is another basis, and another size is \code{sd},
    which a held coefficient makes exact: the forest's contribution then has
    prior standard deviation \code{sd}. A held forest with no basis that states
    no \code{sd} has the size \code{\link{bart}} gives the one forest of an
    ordinary fit under a continuous response, a quarter of the response's
    range, and one unit of the link's error under \code{"probit"} and
    \code{"logistic"}. A numeric column held enters as given, its forest at
    full size wherever the column is 1: \code{basis = z - 0.5} held is the
    forest times -0.5 and 0.5. A held coefficient cannot adapt to an
    \code{sd} that is too small; where the effect may be larger, state
    \code{sd} or let the coefficient be drawn. \code{fixed} alone is taken as
    \code{fixed()}, and \code{student} alone as \code{student()}; no other
    value is, a response family among them. Outside the argument that takes
    them the constructors are \code{dbartsPriors$fixed()} and
    \code{dbartsPriors$student()}. The prior is stated for a forest of a model
    of several, at creation; a model of one forest has no coefficient to give
    a prior, and its size is drawn through the fitting function's leaf prior,
    as \code{leaf.prior = normal(k = chi(1.5, 2))} or
    \code{normal(sd = invchi(1.5, s))}.
    \code{$getLeafPrior(f)$amplitude} reads it back.
  }
  \item{proposal.probs, n.perturb.cuts, treeShift}{
    This forest's tree moves - its proposal mixture, perturb window and
    level-shifting step - as \code{\link{dbartsControl}} describes each. Each
    one stated here governs this forest and each one left out takes the
    control's: \code{$setControl} changes a forest's setting only where the
    forest states none, and one stated here is fixed for the sampler's life.
    A stated \code{proposal.probs} is a whole mixture, filled against the
    package defaults as \code{dbartsControl()}'s is and not edits to the
    control's, so an omitted \code{birth} is 0.5. A forest frozen by its own
    mixture is still grown by \code{n.grow.sweeps} and replaced by
    \code{warm.start}, as under the control's, and \code{treeShift} acts on
    the constant leaves of two or more trees, as the control's does.
    \code{n =} and \code{p =} each abbreviate more than one argument, so
    write \code{n.trees} and \code{power} in full.
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
  which is the shipped single-forest binary default's own 0.239. A forest's
  default does not depend on how many forests the model has, so each forest
  added widens the prior on the combined location. When every forest carries
  a factor the induced prior standard deviation of the index is
  \eqn{1.049\sqrt{K}}{1.049 sqrt(K)} latent units, 1.484 at
  \eqn{K = 2}{K = 2}, which is 0.989 of the classic \eqn{k = 2}{k = 2} binary
  leaf-scale budget; a numeric column at its default contributes 1.484 for
  each standard deviation of the column the row is from 0. When one forest
  carries none, its amplitude is Cauchy and has no variance to enter that
  budget with, and the others contribute
  \eqn{1.049\sqrt{K - 1}}{1.049 sqrt(K - 1)}, 0.699 of the budget at
  \eqn{K = 2}{K = 2}. Where many forests share one index, state their
  \code{sd}. Read the values in force off
  \code{$getLeafPrior(f)}'s \code{sd} entry, in the response's units.

  Both forests' leaf scales come from the model's own calibration map rather
  than from the leaf prior, which is why a \code{k} hyperprior, a non-default
  \code{k}, a named leaf-prior \code{sd}, and a linear or Gaussian-process
  leaf prior are refused when a second forest is declared, the size of the
  forest with no basis being written on the forest, as
  \code{forest(x1 + x2, sd = 2)}, in the response's units as
  \code{normal(sd = )} is. In a model of several forests this \code{sd}
  states one forest's share of the combined location's prior, per unit of a
  numeric basis, where \code{normal(sd = )} states a single forest's whole
  spread. Every value here is validated at fit time, and
  anything today's engine cannot honor is refused there by name rather than
  dropped.

  In a \code{\link{bart}} or \code{\link{dbarts}} formula each
  \code{forest()} term is one forest of the model, without a separate
  \code{forests =} list:
  \code{y ~ forest(x1 + x2) + forest(x1 + x2, basis = z)}. A term stands at
  the top of the formula's right-hand side, joined to the others by \code{+};
  a \code{forest()} crossed with another term, as \code{z:forest(x1 + x2)} or
  \code{z * forest(x1 + x2)}, is refused with the \code{forest()} to write in
  its place. The forest with no \code{basis} is the forest with no multiplier
  and the model's first, wherever it is written; a forest's default tree
  count and tree prior go by whether it has a \code{basis} and not by where
  it is written. It may be left as plain
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
# its coefficients' prior variance is 1/2
indexSd <- sqrt(calibration$k.scale^2 * 0.5 * rowSums(basis^2))

# the same model written in a formula, the first forest written out; the
# columns of a basis are named as lm() would name them
d <- data.frame(y = y, x1 = x[, 1], x2 = x[, 2], x3 = x[, 3], z = z,
                w = 50 + 10 * rnorm(n))
written <- dbarts(y ~ forest(x1 + x2 + x3) +
                    forest(x1 + x3, basis = factor(z), n.trees = 25L, sd = 1.5),
                  d, control = dbartsControl(n.chains = 1L, n.trees = 25L))
colnames(written$data@bases[[2L]])

# a standardized multiplier: predict centers and scales the new rows by the
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

# two columns are a forest each, with a coefficient each; I() holds arithmetic
perColumn <- dbarts(y ~ forest(x1 + x3) + forest(x1, basis = z) +
                      forest(x1, basis = I(w / 50)),
                    d, control = dbartsControl(n.chains = 1L, n.trees = 10L))
lapply(perColumn$data@bases[-1L], colnames)
}
\keyword{models}
