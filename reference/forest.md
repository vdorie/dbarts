# Forest Specification for Multi-Forest Models

One forest of a model of several forests. Such a model fits the mean as
a sum of forests, each multiplied by columns of the data:
`y ~ forest(x1 + x2) + forest(x1 + x2, basis = dose)` is \\a_1 f_1(x) +
a_2 \cdot \mathrm{dose} \cdot f_2(x)\\, where each \\f\\ is a forest and
each \\a\\ a coefficient the sampler draws. A forest is written as a
term of a [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) or
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) formula,
as in that example, or as an element of the `forests` list of
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) or
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md):
`forests = list(forest(), forest(basis = dose))`. Either way `vars` is
the predictors the forest splits on and `basis` is what it is multiplied
by. Every setting is per forest and stated on the forest.

The constructor is not exported: it resolves by bare name inside the
arguments that take it, and elsewhere is written
`dbartsForests$forest(...)`; see
[`dbartsForests`](https://vdorie.github.io/dbarts/reference/dbartsForests.md).

## Usage

``` r
forest(
    vars = NULL, basis = NULL, sd = NULL,
    n.trees = NULL, base = NULL, power = NULL, amplitude = NULL,
    interactions = NULL, blocks = NULL,
    proposal.probs = NULL, n.perturb.cuts = NULL, treeShift = NULL)
```

## Arguments

- vars:

  The predictors this forest splits on. In short:

  - It is the one argument given without its name: `forest(x1 + x2)`.
    Every other argument is given by name.

  - Left out, the forest splits on every predictor of the fit.

  - In a formula it is terms joined by `+`, as on the right of any model
    formula: `forest(x1 + x2)`, `forest(log(x1) + factor(g))`,
    `forest(. - z)`.

  - In a `forests` list the same terms select among the fit's
    predictors, `forest(x1 + x3)`, `forest(. - x2)`; names or positions
    do too, `forest(c("x1", "x3"))`, `forest(c(1, 3))`.

  - A name of a predictor is that predictor. Anything else is read where
    `forest()` is called, at that moment.

  The rest of this item is the detail of each.

  In a formula given to
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) or
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) it is
  written as the right-hand side of a model formula of the forest's own:
  `forest(x1 + x2)`, `forest(log(x1) + factor(g))`, `forest(. - z)`,
  where `.` is every column of `data` but the response and `-` removes a
  term. The fit's predictors are the terms of all of its forests, each
  once, and a forest splits on the columns of its own terms; a term that
  gives several columns, such as `poly(x1, 2)` or a factor under
  `factors = "indicators"`, is all of them. Predictors may also be
  named, as `forest(c("x1", "x2"))`: a name of a column of `data` is the
  term that name written out is, and any other name is a column of the
  predictors as the fit builds them, one of a factor's indicator columns
  for instance. A position is not taken in a formula. An
  [`offset()`](https://rdrr.io/r/stats/offset.html) and an intercept
  term (`1`, `0` or `- 1`) belong to the fit and are written beside the
  forests, not inside one. A removal written beside the forests takes
  its term from the forest with no `basis` and from no forest that names
  predictors of its own: the forest with no `basis` has the terms of the
  formula with its `forest()` replaced by what is inside it and every
  forest with a `basis` left out. So
  `y ~ forest(.) - x3 + forest(x3, basis = z)` and
  `y ~ . - x3 + forest(x3, basis = z)` are one model, in which the
  second forest alone splits on `x3`; a removal that leaves the first
  forest no predictor is refused. Written with no first argument, that
  forest is every predictor the other forests name, less what is removed
  beside it. A removal of a term the forest with no `basis` does not
  have is passed by, as R passes it by in a plain formula:
  `y ~ x1 + x2 + forest(x3, basis = z) - x3` is the model without the
  removal. [`update`](https://rdrr.io/r/stats/update.html) cannot take a
  term out of a `forest()`: `update(f, . ~ . - x3)` returns `f`
  unchanged when `x3` stands inside one.

  In a `forests` list it selects among the fit's predictors, in one of
  two ways. When a name in it is a predictor of the fit, it is read as
  terms in the same way: `forest(x1 + x3)`, `forest(. - x2)`,
  `forest(. - log(x1))`. There `.` is every column of the fit's
  predictors, and a term names a predictor as the fit holds it, by its
  term label when the fit has a formula and by its column name when it
  has a matrix, written as code or as a backticked name; a term that
  names no predictor is refused, a removed one too. Anything else is a
  value: it is evaluated once, where `forest()` is called and at that
  moment, and gives the predictors' names or their positions among the
  fit's columns: `forest(c("x1", "x3"))`, `forest(c(1, 3))`, or a
  variable holding either. A warning that evaluation raises is kept and
  raised when a fit uses the value. A name or a position given twice is
  refused, and so are a logical, a factor and a formula: the multiplier
  of a forest is its `basis`, by name. A column with no name is selected
  by position or by `.` only, and a name that two columns share is
  refused. A number is a position only as the whole argument: among
  terms, as in `forest(. - 2)`, it is refused, and a column named `2` is
  written backticked. With a factor expanded to indicator columns,
  removing one of them, as `forest(. - g.u)`, leaves the forest able to
  separate that level through the others; remove the factor,
  `forest(. - g)`, to keep it out.

  A name of a predictor is that predictor, and anything else is read
  where `forest()` is called, at that moment. A forest built in a loop,
  by `lapply(names, forest, basis = z)` or by
  [`Map()`](https://rdrr.io/r/base/funprog.html) therefore keeps the
  selection each call was given, whatever its variables hold by the time
  the model is fitted, and a forest saved to a file carries its
  selection with it. A variable that does not exist yet when `forest()`
  is called is not looked for later. One case needs care: a predictor
  hides a variable of yours with the same name, so in `forest(x1)` with
  a predictor `x1` and a variable `x1` that holds names, the predictor
  is meant. Where the two could meet, hand the value over:
  `do.call(forest, list(vars = nms))`.

  Any forest may be restricted, a single declared forest included: given
  the same residual scale estimate, that fit is the fit on the named
  columns alone, and its other columns report no splits. The estimate is
  taken from every column unless `sigest` states it, so the two fits
  agree draw for draw when it does. On a single forest a `dart` tree
  prior beside `vars` lays its Dirichlet over the named columns, every
  other column reporting probability 0, and `split.probs` keeps its
  ratios among the named columns and must give one of them a positive
  probability (see
  [`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)).
  On every forest of a hazard fit the `period` column the fit appends is
  always allowed, named or not: `vars` restricts the columns the caller
  supplied. The restriction holds for the sampler's life; `$setModel`
  refuses a model that states another, or none.

- basis:

  What the forest is multiplied by. In short:

  - A name of a column of `data` is that column: `basis = dose`.

  - A basis is one numeric column, or a factor, a character column or a
    logical column of two levels, with a coefficient for each level. A
    number is one column, whatever its values: `basis = z` for a `z` of
    0s and 1s is one coefficient, and `basis = factor(z)` is the two
    levels.

  - `I(dose + age)` is one column, their sum. Several columns,
    `dose + age`, and three or more levels are refused: give each
    column, or each level but the first, a forest of its own.

  - Anything that is not a column of `data` is read where `forest()` is
    called, at that moment, with or without a tilde written in place:
    `forest(basis = w)`, `forest(basis = ~ w)`.

  - A formula held in a variable, `b <- ~ w` and then
    `forest(basis = b)`, is read the same way, when `forest()` is
    called, and is left as you made it.

  - [`cbind()`](https://rdrr.io/r/base/cbind.html) is refused: write its
    one column alone, or give each column a forest of its own.

  The rest of this item is the detail of each.

  **What it does.** With a basis of columns \\B_1, \ldots, B_K\\ the
  forest \\f\\ enters the fit as \\(\sum_k a_k B_k(x_i)) f(x_i)\\, each
  column with a coefficient \\a_k\\ of its own, called its amplitude:
  \\K = 1\\ for a numeric column, and the two columns of a factor's
  levels, \\K = 2\\. A forest with no basis is multiplied by a single
  coefficient. A model has one forest with no basis, or none; every
  other forest states one, and a model whose only forest has a basis is
  refused.

  **Writing it.** A basis is the right-hand side of a model formula,
  without the tilde, and means what it means in
  [`lm`](https://rdrr.io/r/stats/lm.html):

  - a name is a column: `basis = dose`;

  - `+` separates columns: `dose + age` is two columns, which this
    version refuses (below);

  - arithmetic is written inside
    [`I()`](https://rdrr.io/r/base/AsIs.html): `I(dose + age)` is one
    column, the sum, and `I(dose / 30)` is `dose` rescaled;

  - a function of columns is a term, as in `lm`: `log(dose)`,
    `scale(age)` and `poly(dose, 1)`;

  - a factor, a character column or a logical column of two levels gives
    one column for each, none left out: `factor(z)` for a `z` of 0 and 1
    is the two columns whose coefficients are \\(b_0, b_1)\\, in the
    order of the levels;

  - `dose:age` is one column, the product;

  - there is no constant column unless `1` asks for one, and `1 + dose`,
    a column of ones beside `dose`, is refused, while `0 + dose` and
    `dose - 1` are `dose`.

  A basis is either one factor (or character or logical column) or
  numeric columns; the two are not mixed in one forest. The forest does
  not center or scale a basis, so a multiplier that should be
  standardized is written so: `basis = scale(age)`. A numeric column
  constant over the observations in the likelihood has no standard
  deviation to take the default `sd` from, and is refused unless `sd` is
  stated.

  **One size.** In this version a basis is one numeric column, or a
  factor, a character or a logical vector of two levels, counted on the
  rows the fit keeps once a level they leave empty is dropped. Three or
  more levels are refused, with the form to write, a forest for each
  level but the first:
  `forest(x1, basis = arm == "b") + forest(x1, basis = arm == "c")`.
  Several numeric columns are refused however they are written,
  `dose + age`, `poly(dose, 2)`, `splines::ns(dose, 3)`, a matrix, or a
  factor's indicator columns handed over as numbers, with the form to
  write, a forest for each column:
  `forest(x1, basis = dose) + forest(x1, basis = age)`; indicator
  columns are handed over as the factor itself. `1 + dose` is refused
  for its constant column: write `basis = dose`. Each is refused drawn
  or held, at every place a basis is given, and a basis handed over as a
  value is told the form in words.

  **Which kind it is.** The class of the basis decides whether the
  forest multiplies levels or numbers, the same way at every place a
  basis is given: `forest()` written as code or with a value,
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)'s
  `bases`, `$setForestBasis` and `predict`'s `bases`. A factor, ordered
  or not, a character vector and a logical vector are levels, with a
  coefficient for each: in the order of a factor's levels, in sorted
  order for a character vector, `FALSE` then `TRUE` for a logical
  vector. An integer or a double is numbers, whatever its values, and a
  numeric matrix is its columns. The values are never looked at to
  decide: a column of 0s and 1s is one column with one coefficient, and
  is two levels only when written `factor(z)` or as a logical vector. A
  matrix of characters or of logicals, and a factor matrix, is refused,
  one of a single column included, since a vector is levels and a matrix
  is numbers: hand over its column as a vector, `cm[, 1]`. An array of
  one dimension is its vector. A data frame and any other class are
  refused too. Handed over as a value, the refusal names the class;
  written as code, a class a model frame does not take, a data frame or
  a list, is refused in R's own words.

  The kind is fixed when the forest is created. A forest created with no
  basis cannot be given one, and a forest of levels keeps the same
  levels, in the same order, for as long as the sampler lives, so that
  coefficient \\j\\ multiplies level \\j\\'s rows after any replacement.
  See `$setForestBasis` in
  [`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  and the `bases` argument of `predict` in
  [`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md) for
  what a replacement takes.

  A level that no row the fit keeps has gives no column and no message,
  written as code or handed over as a value: `factor(z, levels = 0:2)`
  of a `z` with no 2, or a factor of three levels under a `subset` that
  leaves two of them, is the two columns of the levels that remain. The
  count is made after the rows are cut, and fewer than two levels left
  is refused. A row of weight zero is a row the fit keeps.

  A missing value is looked for on the rows the fit keeps, after
  `subset` and the `na.action`, as `lm` looks for one in its variables.
  On a kept row it is refused. On a row the fit drops it is not looked
  at, so a basis fits the same way written as code or handed over as a
  value. `NA` is not a level: a factor made with
  [`addNA()`](https://rdrr.io/r/base/factor.html) or `exclude = NULL`
  that has a kept row at `NA` is refused, and the remedy is to give that
  level a name, as in `levels(f)[is.na(levels(f))] <- "none"`. With no
  kept row at it the `NA` level is dropped like any other level no row
  has.

  Some forms a model formula would take are refused, each by name and
  with what to write instead: `*` between terms (write `I(dose * age)`
  for the product, or give `dose`, `age` and `I(dose * age)` a forest
  each for all three); a term multiplied or divided by a number (write
  `I(dose / 30)`); `-` between terms, other than `- 1`; `^`, `/`, `%in%`
  and `|` between terms; `.`; an
  [`offset()`](https://rdrr.io/r/stats/offset.html);
  [`cbind()`](https://rdrr.io/r/base/cbind.html) (write its one column
  alone, or give each column a forest of its own); and a term that calls
  one of `normal`, `fixed`, `student`, `cauchy`, `linear`, `gp`, `cgm`,
  `dart`, `chisq`, `chi`, `invchi`, `forest` or `varianceForest`, since
  a prior is not stated on a term of a basis.

  **Names of the columns.** They are named as
  `coef(lm(y ~ 0 + <basis>))` names them: `dose`; `factor(z)0`,
  `factor(z)1`; and, for a basis of several columns, which is read and
  named before it is refused and whose refusal quotes the names, `dose`,
  `age`, `poly(dose, 2)1`, `poly(dose, 2)2` and `(Intercept)` for the
  column of ones. They are the column names of `data@bases[[f]]` on the
  sampler. A forest's label, where its place in a `forests` list has no
  name, is the text of its basis, `"scale(age)"`, or `forest<i>` for a
  forest with no basis or with one given as a value; see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md).
  Wherever a method or an `extract` or `predict` takes a `forest`, a
  string is a forest's label (failing an exact match, the label that is
  the same code, so `"I(dose / 30)"` finds `I(dose/30)`), a number is
  its position as base R reads a list element by number, and `forest<i>`
  names position *i*, though a per-forest margin of a fit is named by
  the forests' labels. A string is never a position: `"2"` is a forest
  labeled `"2"` or none. A string that is one forest's label and another
  position's name, or the same code as two labels, is refused.

  **Rows.** A basis has a value for every row of `data`, and is cut to
  the rows the fit keeps: `data` and `subset` are each evaluated once,
  for the fit, and the basis is given the same rows in the same order,
  less those the `na.action` drops. So with
  `data = d[sample(nrow(d)), ]` or `subset = sample(n, 100)` the basis
  has the rows the fit has.
  [`scale()`](https://rdrr.io/r/base/scale.html),
  [`poly()`](https://rdrr.io/r/stats/poly.html) and anything else a term
  computes across rows, `I(age - mean(age))` among it, are computed on
  every row of `data`, as in `lm`, whatever `subset` keeps. A level of a
  factor that no kept row has gives no column. Refused are a missing
  value in a row the fit keeps, where `lm` would drop the row; a factor,
  a character column or a logical column with one level on the rows
  kept; and a numeric column that is zero on all of them.

  **New rows.** `predict` builds the basis again from the new rows, as
  `predict` does for `lm`. Taken from the fitted rows are the center,
  scale and knots of [`scale()`](https://rdrr.io/r/base/scale.html),
  [`poly()`](https://rdrr.io/r/stats/poly.html), `ns()` and `bs()`, and
  the levels of a factor. Everything else is computed again on the rows
  given: `I(age - mean(age))` uses the mean of the new rows, and
  `factor(dose > median(dose))` their median; see
  [`SafePrediction`](https://rdrr.io/r/stats/makepredictcall.html).
  Every column of `data` the basis names must be among the new rows, a
  factor's columns keep the order of the fit's levels, and a level the
  fit did not have is refused. A basis given as a value cannot be built
  again: give it at the new rows through `predict`'s `bases` argument.

  **Where a name is found.** It depends on who wrote the code.

  In a call of `forest()` that you write, as in a `forests` list or for
  a forest built ahead of the fit, a name of a column of `data` is that
  column. Every other variable of yours that the code uses is copied
  when `forest()` is called. This is so whether or not a tilde is
  written in place: `forest(basis = ~ I(dose + age))` is
  `forest(basis = I(dose + age))`. A forest built in a loop therefore
  keeps what each call was given,


          fs <- list(forest())
          for (k in c(10, 30))
            fs[[length(fs) + 1]] <- forest(x1, basis = I(dose / k))
          

  is one forest multiplied by `dose / 10` and one by `dose / 30`, at the
  fit and at `predict`, whatever `k` holds later, and a saved forest
  carries its numbers with it. A variable that does not exist when
  `forest()` is called is not looked for later, and is refused. What R
  or an attached package supplies, a function such as `scale` or a
  constant such as `pi`, is not copied: a basis that names a column of
  `data` looks it up when the model is fitted. With no data frame, as
  with the matrix interface and
  [`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md),
  there is no column to name and every variable is yours.

  One case needs care: a column of `data` hides a variable of yours with
  the same name, whether that variable holds a column, a set of names or
  a one-sided formula. Where the two could meet, hand the value over:
  `do.call(forest, list(basis = w))` for a column you hold,
  `list(basis = as.name(nm))` for a column of `data` by its name,
  `list(basis = f)` for a one-sided formula, and `list(vars = nms)` for
  predictors by name.

  A formula made elsewhere, held in a variable or handed over as above,
  is read in the same way and at the same moment, when `forest()` is
  called: a name of a column of `data` is that column, and every other
  variable the formula uses is copied, as it is then, from where the
  formula was made. What counts is when `forest()` is called, not when
  the formula was made. In


          fs <- list(forest())
          for (k in c(10, 30)) {
            f <- ~ I(dose / k)
            fs[[length(fs) + 1]] <- forest(x1, basis = f)
          }
          

  each forest has its own `k`, because `forest()` is called inside the
  loop. Formulas made in a loop and given to `forest()` only after it
  all see the last `k`: a formula holds its environment and no value, as
  everywhere in R. A variable changed after `forest()` is called changes
  neither the fit nor a prediction. The formula itself is not touched:
  it is the formula you made, in its environment, and nothing is
  assigned there.

  A `forest()` term inside the formula of
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) or
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) is
  part of that formula and is read with it, as R reads a formula, with a
  tilde or without: in `data` and then where the formula was written,
  when the model is fitted and again by `predict`. Change a number such
  a term uses between the fit and the prediction and the prediction
  changes, with no message.

  An argument that [`lapply()`](https://rdrr.io/r/base/lapply.html),
  [`Map()`](https://rdrr.io/r/base/funprog.html) or another function of
  base R writes for you, as in `lapply(columns, forest, vars = "x1")`,
  is a value handed over, whatever the data's columns are called. A name
  or a call that such a function is handed and passes on, as in
  `Map(forest, nms, MoreArgs = list(basis = quote(dose)))`, is code, as
  it is through [`do.call()`](https://rdrr.io/r/base/do.call.html).

  **A value.** What is handed over as an object, through
  [`do.call()`](https://rdrr.io/r/base/do.call.html) or in a call built
  by a program, is used as it is: a numeric vector or matrix is its
  columns, and a factor, a character vector or a logical vector is its
  levels, as written as code. Its column keeps the name it has. Its
  forest's label is `forest<i>`. It has one row for every row of `data`
  and is cut to the rows the fit keeps; a value already cut to those
  rows is refused by name. A level that no row the fit keeps has gives
  no column, whether `data`, `subset` or the `na.action` is what empties
  it. A value that the rows kept leave with one level, or one numeric
  column that is zero on all of them, is refused as the same basis
  written as code is. A single string or number is refused:
  `basis = "dose"` names no column, so write `basis = dose`.

  A basis is stated here or on the data object, not both:
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) given
  an already-built `dbartsData` refuses this argument, the data object's
  own `bases` being where it is stated then, while
  [`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md)
  takes a data object and installs the bases declared here in place of
  those it carried. After creation a basis is replaced by
  `$setForestBasis`, whose formula is read at once, where it was made;
  see
  [`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md).

  A lone forest carrying a `basis` is refused, the coefficients being
  what distinguish a forest from another. For varying coefficients
  declare an intercept forest plus one basis forest per covariate,
  `forests = list(forest(), forest(basis = z1), forest(basis = z2))`, or
  use one forest with `linear()` leaves.

- n.trees, base, power:

  This forest's tree count and tree-structure prior. Each one stated
  here governs this forest, and each one left out takes the default of
  the forest's kind, whatever the forest's place. A forest with a
  `basis` has 50 trees at `base = 0.25`, `power = 3` - shallower and
  fewer than a prognostic forest's, the modulating surface normally
  being the smoother of the two. The forest with no `basis` has the
  fitting function's:
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s
  `n.trees` or the control's, and `tree.prior`. A value stated here on
  that forest governs it over the control's `n.trees` and the tree
  prior's `base` or `power`, being the more specific of the two
  declarations, while
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) refuses
  its own `n.trees` beside one stated here, the two being the same
  count. Where every forest has a `basis` the fitting function's tree
  count and tree prior have no forest to belong to, and one that is
  stated is refused, the default's own value included: state it on a
  forest. A sampler's `control@n.trees` and `model@tree.prior` hold its
  first forest's.

- sd:

  The size of what this forest contributes, in the response's units: the
  units of \\y\\ for a continuous response, and of the latent index
  under `"probit"` and `"logistic"`, where the link's own error has
  standard deviation 1 and \\\pi/\sqrt{3}\\. One number; a name on it is
  ignored. What it is the size of depends on what multiplies the forest.

  **No basis.** The forest enters as \\a\\f(x)\\, and `sd` is the size
  of \\a\\f(x)\\.

  **A factor** (a factor, a character vector or a logical vector, of two
  levels). The forest enters with one coefficient for each level, and
  `sd` is the size of the difference between the two: the treatment
  effect where the levels are control and treated.

  **A number** \\w\\. A column that is an ordinary covariate should be
  centered: `basis = scale(w)`, or `basis = I(w - c)` for a reference
  value \\c\\. The other forests describe the response where \\w\\ is 0,
  and a column that never comes near 0, an age or a year, makes them an
  extrapolation and the fit poor, with nothing in the chain to show it.
  A dose or an exposure, where 0 means none, is left as it is. The
  forest enters as \\w\\(a\\f(x))\\, and `sd` is the size of
  \\a\\f(x)\\, the change in the response for one unit of \\w\\; the
  term at an observation has size \\\|w\|\\ `sd`. The column is used as
  given: it is never centered and never divided, so a stated `sd` is per
  unit of \\w\\, and the same `sd` with \\w\\ in grams and in milligrams
  states two priors a thousand times apart.

  Where the coefficient is held (`amplitude = fixed()`) it is 1, what it
  multiplies is the basis, and "size" is the prior standard deviation,
  exactly. Where it is drawn the size is itself uncertain, and `sd` is
  its prior median.

  Not stated, `sd` for a forest with no basis is twice
  \\\mathrm{sd}(y)\\ under a continuous response where its coefficient
  is drawn, the size
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) gives its
  forest (a quarter of the response's range) where it is held, and one
  unit of the link's error under `"probit"` and `"logistic"` either way;
  \\\mathrm{sd}(y)\\, or that unit, for a factor; and for a number that
  same size for one standard deviation of the column, \\\mathrm{sd}(y) /
  \mathrm{sd}(w)\\ per unit, held or drawn. None of these depends on how
  many forests the model has. A default fit is then the same fit
  whatever units \\w\\ is in. \\\mathrm{sd}(y)\\ and \\\mathrm{sd}(w)\\
  are taken once, over the observations in the likelihood when the
  sampler is created (those `subset` and the missing-value handling
  keep, with a positive weight, and not masked out), and kept: replacing
  the column, the weights or the mask later does not move the prior,
  unless the column is replaced with
  `$setForestBasis(updateBasisScale = TRUE)`, and a copy, a saved and
  reloaded sampler and a sampler made again from its own control, model
  and data keep it. Where those observations hold fewer than two
  distinct values of the response, both are taken over every
  observation, under one warning (two by the route of
  [`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md)
  and `new()`, see its `active`). `"probit"` and `"logistic"` read no
  \\\mathrm{sd}(y)\\, and where fewer than two observations are in,
  \\\mathrm{sd}(w)\\ is taken over every one, under one warning; two or
  more that hold one value of \\w\\ are refused. \\\mathrm{sd}(y)\\ is
  of the response net of its offset; on a gaussian response that does
  not vary over those observations it is a quarter of the window of
  width one the response is then scaled by. \\\mathrm{sd}(w)\\ is kept
  with the data, in its `basis.scale` slot (see
  [`dbartsData`](https://vdorie.github.io/dbarts/reference/dbartsData.md)),
  and `$getLeafPrior(f)$basis.scale` reports it. A caller who replaces a
  column during burn-in and wants the default of the new one replaces it
  with `$setForestBasis(f, w.new, updateBasisScale = TRUE)`, which takes
  \\\mathrm{sd}(w)\\ again over the observations in the likelihood then
  (see
  [`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)).
  A column that is constant over those observations has no standard
  deviation and needs `sd` stated.

  A treatment belongs in a factor or a logical vector, where `sd` is the
  treatment effect. Written as a 0/1 number it is scaled as any number:
  \\\mathrm{sd}(z)\\ is 0.5, 0.4 and 0.22 at 50, 20 and 5 percent
  treated, so the default effect is 2, 2.5 and 4.6 times the one
  `factor(z)` gets; `sd = sd(y)` on the number states the factor's. The
  same holds for a coding written as a number and held,
  `basis = z - 0.5` with `amplitude = fixed()`.

  A held coefficient cannot adapt to an `sd` that is too small. At the
  default it is sized for a slope of about one \\\mathrm{sd}(y)\\ per
  \\\mathrm{sd}(w)\\; where the effect may be larger, under `"probit"`
  and `"logistic"` above all, state `sd` or let the coefficient be
  drawn.

  A default is not a stated number: `$getLeafPrior(f)$sd` reports the
  number of the response's units it comes to, a numeric forest's named
  by its column, and `$getLeafPrior(f)$sd.stated` is `FALSE`. A stated
  number is kept as written and divided by the unit once, at each
  creation. A control or a `forests` list taken to other data states the
  same numbers in that data's units, and the sampler made from it
  derives its own. It is one number, positive and finite. A vector of
  any other length and anything that is not a number are refused. The
  one forest of a model of one forest takes it as the fitting function's
  `leaf.prior = normal(sd = )`, in the same units, and stating the size
  at both places is refused.

  It is the one argument a live sampler restates:
  `$setLeafPrior(forests = list(forest(sd = ), ...))` states it from
  then on, also where it is the number in force, and leaves the recorded
  \\\mathrm{sd}(w)\\. `$getLeafPrior(f)$leaf.prior` returns the
  statement, `forest(sd = )` where one was stated and `forest()` where
  none was, and `forest()` written back leaves its forest as it is (see
  [`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)).
  A forest's kind is fixed at creation: `$setForestBasis` refuses to
  give a basis to a forest created without one, and moves no forest's
  prior.

  `print` of a fit and `show` of a sampler give a short block of one
  line for each forest of a model of several: its label, what it
  multiplies, for a numeric column the mean and sd it has over the
  observations in the likelihood as they stand, and the sd the forest
  resolved to. For a numeric column that is its sd per sd of the column
  as it now stands, named by the forest's label, the size of the
  forest's term over one sd of the column in the response's units: for a
  default it is \\\mathrm{sd}(y)\\ (the link's error sd under `"probit"`
  and `"logistic"`) wherever the column, and the observations it is read
  over, are those the default was taken from, and for a stated `sd` it
  is that `sd` times the column's sd. Where the column is constant it is
  the `sd` per unit. For a factor it is the sd between the two levels.
  The sd per unit of a numeric column, which is what `forest(sd = )`
  states, is in the full block and in the reader. The last line of the
  short block says where each forest's prior is to be had in full:
  [`summary()`](https://rdrr.io/r/base/summary.html) of a fit, which
  prints the full block, or `$getLeafPrior()` of a sampler, which
  returns it. The `summary` of a fit and the creation of a sampler with
  `verbose = TRUE` give the full block. Its header gives the unit the
  sds are in: \\\mathrm{sd}(y)\\, as taken above, or the link's error
  sd. A numeric line gives the sd per unit and per sd of the column and
  adds the column's range; every line adds whether the sd was stated or
  is the default, and whether the coefficient is held. A default whose
  column has another sd over the observations now in the likelihood than
  the one it was taken from (the column replaced, the weights or the
  active rows changed, or, under `"gaussian"`, every observation read
  because those in the likelihood held fewer than two values of the
  response) says so. A fit made with `keepFits = FALSE` keeps no bases
  or labels and gives no block. Two things to look for. The mean of a
  numeric column, in either block, and its range, in the full one: a
  mean far from 0 is a column that never comes near 0, an age or a year,
  and should be centered, unless 0 means none, as for a dose; and one
  wild value shows in the range, and moves the mean and sd. And the
  figure per sd of the column, in either block, the size of the forest's
  term over one sd of the column (`sd` times the column's sd): for a
  default it is \\\mathrm{sd}(y)\\ on every numeric line whose column,
  and the observations it is read over, are those the default was taken
  from, so a line that differs has had its column or its observations
  changed, and for a stated `sd` it is where a number in the wrong units
  shows. Nothing in the fit warns of either.

- interactions, blocks:

  Optional
  [`interactions`](https://vdorie.github.io/dbarts/reference/interactions.md)
  and [`blocks`](https://vdorie.github.io/dbarts/reference/blocks.md)
  constraints on this forest. The arguments of the same names on
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) are
  those of the forest with no `basis`, wherever it stands, so one stated
  there and on that forest is refused as one constraint given twice. A
  forest with a `basis` states its own here, and where every forest has
  one those arguments of the fitting function are refused. Holding the
  two forests to different structures is the calibrated-additivity
  idiom - an additive or low-order modulating forest beside a free
  prognostic one. A `blocks` partition covers the columns the forest may
  split on, i.e. the `vars` subset when one is given, and the forest's
  own tree count.

- amplitude:

  The prior on this forest's coefficient: the one a forest with no
  `basis` carries, or one for each level of its `basis`. Left out, it is
  drawn each sweep. `fixed()` holds it for the sampler's life, the
  choice being made at creation, wherever the forest stands: a held
  coefficient is 1, and what it multiplies is the basis, so a forest
  with no `basis` enters at 1, a forest on a numeric column at the
  column as given, and a forest on a factor, a character or a logical
  vector of two levels at 0 for the first level and 1 for the second.
  The forest is then held out of the first level and enters at full size
  for the second, so it is the second level's difference from the first.
  Another coding is another basis, and another size is `sd`, which a
  held coefficient makes exact: the forest's contribution then has prior
  standard deviation `sd`. A held forest with no basis that states no
  `sd` has the size
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) gives the
  one forest of an ordinary fit under a continuous response, a quarter
  of the response's range, and one unit of the link's error under
  `"probit"` and `"logistic"`. A numeric column held enters as given,
  its forest at full size wherever the column is 1: `basis = z - 0.5`
  held is the forest times -0.5 and 0.5. A held coefficient cannot adapt
  to an `sd` that is too small; where the effect may be larger, state
  `sd` or let the coefficient be drawn. `fixed` alone is taken as
  `fixed()`, and no other value is. Outside the argument that takes it
  the constructor is `dbartsPriors$fixed()`. It is stated for a forest
  of a model of several; a model of one forest has no amplitude to hold.

- proposal.probs, n.perturb.cuts, treeShift:

  This forest's tree moves - its proposal mixture, perturb window and
  level-shifting step - as
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  describes each. Each one stated here governs this forest and each one
  left out takes the control's: `$setControl` changes a forest's setting
  only where the forest states none, and one stated here is fixed for
  the sampler's life. A stated `proposal.probs` is a whole mixture,
  filled against the package defaults as
  [`dbartsControl()`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)'s
  is and not edits to the control's, so an omitted `birth` is 0.5. A
  forest frozen by its own mixture is still grown by `n.grow.sweeps` and
  replaced by `warm.start`, as under the control's, and `treeShift` acts
  on the constant leaves of two or more trees, as the control's does.
  `n =` and `p =` each abbreviate more than one argument, so write
  `n.trees` and `power` in full.

## Details

With two forests, the second carrying a two-level factor basis, this is
the Bayesian causal forest \\y = a\\\mu(x) + b_z\\\tau(x) + \epsilon\\:
a prognostic forest \\\mu\\ over every predictor, a modulating forest
\\\tau\\ over the columns `vars` allows, and the amplitudes \\(a, b_0,
b_1)\\ joining them, read back with `$getForestAmplitudes`. Gaussian,
`"probit"` and `"logistic"` responses; under a latent family the
combination is the index rather than the mean, on the link's own fixed
scale, and `"aft"`, `"ordinal"` and `"nbinom"` are refused at creation
naming what each is missing.

What the defaults put on the combined location, since under a latent
family that location IS the index and no drawn \\\sigma\\ stands between
it and the fitted probabilities. At two forests of the shipped shape -
one carrying no basis - a probit model's prior puts \\P(p \< 0.01
\mathrm{~or~} p \> 0.99)\\ at 0.238, which is the shipped single-forest
binary default's own 0.239. A forest's default does not depend on how
many forests the model has, so each forest added widens the prior on the
combined location. When every forest carries a factor the induced prior
standard deviation of the index is \\1.049\sqrt{K}\\ latent units, 1.484
at \\K = 2\\, which is 0.989 of the classic \\k = 2\\ binary leaf-scale
budget; a numeric column at its default contributes 1.484 for each
standard deviation of the column the row is from 0. When one forest
carries none, its amplitude is Cauchy and has no variance to enter that
budget with, and the others contribute \\1.049\sqrt{K - 1}\\, 0.699 of
the budget at \\K = 2\\. Where many forests share one index, state their
`sd`. Read the values in force off `$getLeafPrior(f)`'s `sd` entry, in
the response's units.

Both forests' leaf scales come from the model's own calibration map
rather than from the leaf prior, which is why a `k` hyperprior, a
non-default `k`, a named leaf-prior `sd`, and a linear or
Gaussian-process leaf prior are refused when a second forest is
declared, the size of the forest with no basis being written on the
forest, as `forest(x1 + x2, sd = 2)`, in the response's units as
`normal(sd = )` is. In a model of several forests this `sd` states one
forest's share of the combined location's prior, per unit of a numeric
basis, where `normal(sd = )` states a single forest's whole spread.
Every value here is validated at fit time, and anything today's engine
cannot honor is refused there by name rather than dropped.

In a [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) or
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) formula
each `forest()` term is one forest of the model, without a separate
`forests =` list: `y ~ forest(x1 + x2) + forest(x1 + x2, basis = z)`. A
term stands at the top of the formula's right-hand side, joined to the
others by `+`; a `forest()` crossed with another term, as
`z:forest(x1 + x2)` or `z * forest(x1 + x2)`, is refused with the
`forest()` to write in its place. The forest with no `basis` is the
forest with no multiplier and the model's first, wherever it is written;
a forest's default tree count and tree prior go by whether it has a
`basis` and not by where it is written. It may be left as plain terms,
`y ~ x1 + x2 + forest(x1 + x2, basis = z)`, which is the same model, but
a formula has one such forest: plain predictor terms beside a `forest()`
with no `basis`, and two such terms, are refused. Alone,
`y ~ forest(x1 + x2)` is the single-forest fit `y ~ x1 + x2`, in every
family that takes a formula. The other forests keep the order they are
written in. `forest()` may also be written there as
`dbartsForests$forest()`. A term's predictors and its basis are read as
code; every other argument is evaluated where the formula was written.
The same basis written in a formula's term and in a `forests` list is
the same columns on the same rows, and so the same model.

## Value

A `dbartsForest` specification object, resolved when a sampler is built.

## See also

[`dbartsForests`](https://vdorie.github.io/dbarts/reference/dbartsForests.md),
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md),
[`interactions`](https://vdorie.github.io/dbarts/reference/interactions.md),
[`blocks`](https://vdorie.github.io/dbarts/reference/blocks.md)

[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s ‘Formula
Terms’ section is the full term grammar - the forests of a formula, a
forest's predictors, and the refusal list - of which the ‘Details’
paragraph above is the constructor-side summary.

## Examples

``` r
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
#> [1] "factor(z)0" "factor(z)1"

# a standardized multiplier: predict centers and scales the new rows by the
# training rows' mean and sd of w
fit <- bart(y ~ forest(x1 + x3) +
              forest(x1 + x3, basis = scale(w), n.trees = 10L),
            d, n.trees = 10L, n.samples = 10L, n.burn = 10L,
            n.chains = 1L, n.threads = 1L, keepTrees = TRUE)
#> estimating the starting sigma by a dense linear regression on 100 rows and 2 columns; see 'sigest'
#> family = "auto": continuous response detected, fitting family = "gaussian"; set 'family' to override
#> 
#> Running BART with numeric y
#> 
#> number of trees: 10
#> number of chains: 1, default number of threads 1
#> tree thinning rate: 1
#> Prior:
#>  k prior fixed to 2.000000
#>  degrees of freedom in sigma prior: 3.000000
#>  quantile in sigma prior: 0.900000
#>  scale in sigma prior: 0.008161
#>  power and base for tree prior: 2.000000 0.950000
#>  use quantiles for rule cut points: false
#>  level fiber gibbs step: auto
#>  proposal probabilities: birth/death 0.60, swap 0.00, change 0.40, perturb 0.00, rule_gibbs 0.00; birth 0.50; perturb window 1
#> data:
#>  number of training observations: 100
#>  number of test observations: 0
#>  number of explanatory variables: 2
#>  init sigma: 0.948004, curr sigma: 0.948004
#> 
#> Cutoff rules c in x<=c vs x>c
#> Number of cutoffs: (var: number of possible c):
#> (1: 100) (2: 100) 
#> forests, sd in the response's units (sd(y) 1.139):
#>   1 forest1: no basis; sd 2.278 (default)
#>   2 scale(w): numeric basis, mean 0, sd 1, range -2.666 to 1.923; sd 1.139 per unit, 1.139 per sd (default)
#> Running mcmc loop:
#> total seconds in loop: 0.000437
#> 
#> Tree sizes, last iteration:
#> [1] 2 3 3 2 3 3 2 2 2 2 
#> 
#> Variable Usage, last iteration (var:count):
#> (1: 10) (2: 4) 
#> DONE BART
#> 
predict(fit, d[1:3, ])
#>              1        2        3
#>  [1,] 2.927474 1.455927 1.616158
#>  [2,] 2.760547 1.257952 1.637542
#>  [3,] 2.992327 1.244774 1.607548
#>  [4,] 2.534323 1.379762 1.449334
#>  [5,] 2.417520 1.363156 1.109988
#>  [6,] 1.753917 1.157411 1.286790
#>  [7,] 3.029302 0.936659 1.071813
#>  [8,] 2.460682 1.142985 1.244382
#>  [9,] 3.112798 1.896018 1.416407
#> [10,] 2.176968 1.477079 1.020467

# forests built in a loop: each keeps the number its own call was given
forest <- dbartsForests$forest
forests <- list(forest())
for (k in c(25, 100))
  forests[[length(forests) + 1L]] <- forest(x1, basis = I(w / k))
looped <- dbarts(y ~ x1 + x3, d, forests = forests,
                 control = dbartsControl(n.chains = 1L, n.trees = 10L))
vapply(looped$data@bases[-1L], function(basis) basis[1L, 1L], 0) * c(25, 100)
#> [1] 58.93674 58.93674

# two columns are a forest each, with a coefficient each; I() holds arithmetic
perColumn <- dbarts(y ~ forest(x1 + x3) + forest(x1, basis = z) +
                      forest(x1, basis = I(w / 50)),
                    d, control = dbartsControl(n.chains = 1L, n.trees = 10L))
lapply(perColumn$data@bases[-1L], colnames)
#> [[1]]
#> [1] "z"
#> 
#> [[2]]
#> [1] "I(w/50)"
#> 
```
