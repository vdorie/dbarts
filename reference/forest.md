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
    amplitude.prior.variance = NULL)
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
  On the single forest of a hazard fit the `period` column the fit
  appends is always allowed, named or not: `vars` restricts the columns
  the caller supplied. A hazard fit with several forests takes each
  forest's `vars` as written, so a restricted forest there splits on
  `period` only when `vars` names it. The restriction holds for the
  sampler's life; `$setModel` refuses a model that states another, or
  none.

- basis:

  What the forest is multiplied by. In short:

  - A name of a column of `data` is that column: `basis = dose`.

  - `+` separates columns, each with a coefficient of its own:
    `dose + age`. `I(dose + age)` is one column, their sum.

  - A factor, a character column or a logical column is one column for
    each of its levels.

  - Anything that is not a column of `data` is read where `forest()` is
    called, at that moment, with or without a tilde written in place:
    `forest(basis = w)`, `forest(basis = ~ w)`.

  - A formula held in a variable, `b <- ~ w` and then
    `forest(basis = b)`, is read the same way, when `forest()` is
    called, and is left as you made it.

  - [`cbind()`](https://rdrr.io/r/base/cbind.html) is refused: write the
    columns with `+`.

  The rest of this item is the detail of each.

  **What it does.** With a basis of columns \\B_1, \ldots, B_K\\ the
  forest \\f\\ enters the fit as \\(\sum_k a_k B_k(x_i)) f(x_i)\\, each
  column with a coefficient \\a_k\\ of its own, called its amplitude. A
  forest with no basis is multiplied by a single coefficient. A model
  has one forest with no basis, or none; every other forest states one,
  and a model whose only forest has a basis is refused.

  **Writing it.** A basis is the right-hand side of a model formula,
  without the tilde, and means what it means in
  [`lm`](https://rdrr.io/r/stats/lm.html):

  - a name is a column: `basis = dose`;

  - `+` separates columns: `dose + age` is two columns;

  - arithmetic is written inside
    [`I()`](https://rdrr.io/r/base/AsIs.html): `I(dose + age)` is one
    column, the sum, and `I(dose / 30)` is `dose` rescaled;

  - a function of columns is a term, as in `lm`: `log(dose)`,
    `scale(age)`, and `poly(dose, 2)`, which is two columns;

  - a factor, a character column or a logical column gives one column
    for each of its levels, none left out: `factor(z)` for a `z` of 0
    and 1 is the two columns whose coefficients are \\(b_0, b_1)\\, in
    the order of the levels;

  - `dose:age` is one column, the product;

  - there is no constant column unless `1` asks for one: `1 + dose` is a
    column of ones and `dose`, while `0 + dose` and `dose - 1` are
    `dose`.

  A basis is either one factor (or character or logical column) or
  numeric columns; the two are not mixed in one forest. The forest does
  not centre or scale a basis, so a multiplier that should be
  standardized is written so: `basis = scale(age)`.

  Some forms a model formula would take are refused, each by name and
  with what to write instead: `*` between terms (write `I(dose * age)`
  for the product, or `dose + age + I(dose * age)` for all three
  columns); a term multiplied or divided by a number (write
  `I(dose / 30)`); `-` between terms, other than `- 1`; `^`, `/`, `%in%`
  and `|` between terms; `.`; an
  [`offset()`](https://rdrr.io/r/stats/offset.html);
  [`cbind()`](https://rdrr.io/r/base/cbind.html) (write `+`); and a term
  that calls one of `normal`, `fixed`, `student`, `cauchy`, `linear`,
  `gp`, `cgm`, `dart`, `chisq`, `chi`, `invchi`, `forest` or
  `varianceForest`, since a prior is not stated on a term of a basis.

  **Names of the columns.** They are named as
  `coef(lm(y ~ 0 + <basis>))` names them: `dose`, `age`;
  `poly(dose, 2)1`, `poly(dose, 2)2`; `factor(z)0`, `factor(z)1`;
  `(Intercept)` for the column of ones. They are the column names of
  `data@bases[[f]]` on the sampler. A forest's label, where its place in
  a `forests` list has no name, is the text of its basis,
  `"dose + age"`, or `forest<i>` for a forest with no basis or with one
  given as a value; see
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md).
  Wherever a method or an `extract` or `predict` takes a `forest`, a
  string is a forest's label (failing an exact match, the label that is
  the same code, so `"I(dose / 30)"` finds `I(dose/30)`), a number is
  its position as base R reads a list element by number, and `forest<i>`
  names position *i*, as it does on every per-forest margin. A string is
  never a position: `"2"` is a forest labelled `"2"` or none. A string
  that is one forest's label and another position's name, or the same
  code as two labels, is refused.

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
  `predict` does for `lm`. Taken from the fitted rows are the centre,
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
  written in place: `forest(basis = ~ dose + age)` is
  `forest(basis = dose + age)`. A forest built in a loop therefore keeps
  what each call was given,


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
  levels. Its columns keep their names when every column has one, and
  have none when any column lacks one; two columns of one name are
  refused. Its forest's label is `forest<i>`. It has one row for every
  row of `data` and is cut to the rows the fit keeps; a value already
  cut to those rows is refused by name. A level that no row of `data`
  has is refused, while one that only `subset` empties keeps a column of
  zeros, a value not being built again on the rows kept. A value that
  the rows kept leave with one level, or one numeric column that is zero
  on all of them, is refused as the same basis written as code is. A
  single string or number is refused: `basis = "dose"` names no column,
  so write `basis = dose`.

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

  This forest's prior scale, in units of the response family's own
  latent scale and per unit of basis row norm: the total \\a\\f(x)\\ or
  \\(b_1 - b_0)f(x)\\ is placed at `sd` of them. The unit is
  \\\mathrm{sd}(y)\\ for a gaussian response, `1` for `"probit"` and
  \\\pi/\sqrt{3}\\ for `"logistic"` - the standard deviation of the
  link's own error law, a latent model having no response standard
  deviation to name. A forest whose `basis` rows have median non-zero
  norm \\c\\ contributes the scale named here, the calibration map
  dividing \\c\\ out, so rescaling a basis column does not silently
  rescale the prior. Which of the two channels carries it depends on
  whether the forest has a `basis`: without one it is the half-Cauchy
  median of the forest's scalar amplitude, with one it scales the
  forest's own leaf prior. The two channels default differently, and
  neither default is a bare constant. A forest with NO basis takes `2`
  under a gaussian response, where the unit is the response's own
  \\\mathrm{sd}(y)\\ and a drawn \\\sigma\\ absorbs the difference, and
  `1` under `"probit"` and `"logistic"`, where the unit is the link's
  fixed error scale and nothing does. A forest WITH a basis takes
  \\\sqrt{2/K}\\ in a model of \\K\\ forests, so that declaring more of
  them does not widen the prior on the combined location without bound.
  \\K = 2\\ is the fixed point of both statements, \\\sqrt{2/2} = 1\\. A
  value declared here overrides its default and keeps its per-forest
  reading at every \\K\\, so `sd = 1` on each basis forest recovers the
  pre-\\K\\-aware model exactly. It is one unnamed number, positive and
  finite. For a basis of several columns that one number is the size
  stated for each of them, so `forest(basis = dose + age, sd = 2)`
  states 2 for the column `dose` and 2 for the column `age`; to size two
  columns differently, rescale one of them in the basis, as
  `I(dose / 30)`. A vector of any other length, a named number and
  anything that is not a number are refused. A model of one forest
  states its size as the fitting function's
  `leaf.prior = normal(sd = )`.

  It is the one argument a live sampler restates:
  `$setLeafPrior(forests = list(forest(sd = ), ...))` writes it in the
  channel creation gave the forest, and `$getLeafPrior(f)$leaf.prior`
  reads it back as `forest(sd = )` (see
  [`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)).
  The channel is fixed at creation: a forest created without a basis
  keeps its half-Cauchy median after `$setForestBasis` gives it one, so
  its `forest(sd = )` round-trips through `$setLeafPrior`, while a fresh
  [`dbarts()`](https://vdorie.github.io/dbarts/reference/dbarts.md)
  given the same bases would read the same `sd` as a leaf-scale factor
  and build a different prior.

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

- amplitude.prior.variance:

  Prior variance of the \\N(0, \cdot)\\ amplitudes on this forest's
  basis columns, default `0.5`. Legal only on a forest given a `basis`:
  a forest without one carries a plain scalar amplitude under the
  engine's half-Cauchy scale-mixture prior, whose median is `sd` rather
  than a variance. It is a free multiplier on the induced prior: the
  prior standard deviation of the combined location at row \\i\\ is
  \\\sqrt{\sum_f s_f^2 v_f \\B_f(i,\cdot)\\^2}\\ over the forests
  carrying a basis, with \\s_f\\ read from `$getLeafPrior(f)$k.scale`
  and \\v_f\\ this argument; a basis-free forest's own term is Cauchy
  and has no standard deviation. The budget that sum sits in is set by
  `sd`, whose default already divides it among the \\K\\ forests, so
  raising this argument raises the total rather than redistributing it.
  Under `"probit"` and `"logistic"` that location IS the latent index
  and \\\sigma\\ is pinned, so nothing in the sampler absorbs a
  mis-scaled basis; under a gaussian response it is in
  \\\mathrm{sd}(y)\\ units and a drawn \\\sigma\\ partly does. Every
  input to that expression is readable off the fitted sampler: \\v_f\\
  is the `amplitude.prior.variance` entry of `$getLeafPrior(f)` and
  \\B_f\\ is `data@bases[[f]]`, so the induced prior can be checked
  against what is in force rather than against what the call asked for.
  See the example below. A lone forest carrying a `basis` is refused,
  the amplitudes being what distinguish a forest from another. For
  varying coefficients declare an intercept forest plus one basis forest
  per covariate,
  `forests = list(forest(), forest(basis = z1), forest(basis = z2))`, or
  use one forest with `linear()` leaves.

- amplitude:

  The law of this forest's amplitudes: the one a forest with no `basis`
  carries, or one for each column of its `basis`. Left out, they are
  drawn each sweep. `fixed()` holds them for the sampler's life, the
  choice being made at creation: a forest with no `basis` at 1, and a
  forest on a factor, a character or a logical vector of two levels at 0
  for the first level and 1 for the second. The forest is then held out
  of the first level and enters at full size for the second, so it is
  the second level's difference from the first. For now a held forest is
  taken in two shapes only: a forest with no `basis` at any place but
  the second, and a basis of two columns as the second forest, two
  numeric columns being held as two levels are, the first at 0, so that
  column is not used. Every other shape is refused when held: a forest
  with no `basis` put second in a `forests` list or in a data object's
  `bases`, a basis of two columns at any other place, and a basis of one
  column or of three or more. `$setForestBasis` refuses a held forest a
  replacement of another width. `fixed` alone is taken as `fixed()`, and
  no other value is. Outside the argument that takes it the constructor
  is `dbartsPriors$fixed()`. It is stated for a forest of a model of
  several; a model of one forest has no amplitude to hold.

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
binary default's own 0.239; before the `sd` defaults above it was 0.376.
The \\\sqrt{2/K}\\ factor is what holds that as \\K\\ grows, and it
holds it in two different senses. When EVERY forest carries a basis the
induced prior standard deviation of the index is 1.484 latent units at
every \\K\\ - 0.989 of the classic \\k = 2\\ binary leaf-scale budget -
because the whole location is then a sum of fixed-variance channels.
When one forest carries none, its amplitude is Cauchy and has no
variance to enter that budget with, so the fixed-variance part is
BOUNDED by 1.484 rather than pinned at it, rising from 0.699 of the
budget at \\K = 2\\ toward 0.989 and never reaching it; without the
factor it would instead grow past twice the budget by ten forests. Read
the values in force off `$getLeafPrior(f)`'s `leaf.scale.factor` and
`amplitude.prior.scale` entries.

Both forests' leaf scales come from the model's own calibration map
rather than from the leaf prior, which is why a `k` hyperprior, a
non-default `k`, a named leaf-prior `sd`, and a linear or
Gaussian-process leaf prior are refused when a second forest is
declared. This `sd` is not the leaf prior's: it states this forest's
share of the combined location's prior, per unit of basis row norm,
where `normal(sd = )` states a single forest's whole spread. Every value
here is validated at fit time, and anything today's engine cannot honour
is refused there by name rather than dropped.

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
#> [1] "factor(z)0" "factor(z)1"

# a standardized multiplier: predict centres and scales the new rows by the
# training rows' mean and sd of w
fit <- bart(y ~ forest(x1 + x3) +
              forest(x1 + x3, basis = scale(w), n.trees = 10L),
            d, n.trees = 10L, n.samples = 10L, n.burn = 10L,
            n.chains = 1L, n.threads = 1L, keepTrees = TRUE)
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
#>  level fibre gibbs step: auto
#>  proposal probabilities: birth/death 0.60, swap 0.00, change 0.40, perturb 0.00, rule_gibbs 0.00; birth 0.50
#> data:
#>  number of training observations: 100
#>  number of test observations: 0
#>  number of explanatory variables: 2
#>  init sigma: 0.948004, curr sigma: 0.948004
#> 
#> Cutoff rules c in x<=c vs x>c
#> Number of cutoffs: (var: number of possible c):
#> (1: 100) (2: 100) 
#> Running mcmc loop:
#> total seconds in loop: 0.000664
#> 
#> Tree sizes, last iteration:
#> [1] 3 2 2 3 2 2 3 3 2 2 
#> 
#> Variable Usage, last iteration (var:count):
#> (1: 7) (2: 7) 
#> DONE BART
#> 
predict(fit, d[1:3, ])
#>              1         2         3
#>  [1,] 2.422435 1.4452523 1.7864288
#>  [2,] 2.417635 0.8071004 1.4025173
#>  [3,] 2.917138 0.8665828 1.3632450
#>  [4,] 2.498875 1.0240818 1.3595190
#>  [5,] 2.645444 1.2366439 1.6110356
#>  [6,] 2.261323 1.3595897 1.0811054
#>  [7,] 2.506594 1.5759485 1.5432470
#>  [8,] 2.491697 1.3255230 1.6222154
#>  [9,] 2.185772 1.4219433 1.4804516
#> [10,] 2.480931 1.3449315 0.8162834

# forests built in a loop: each keeps the number its own call was given
forest <- dbartsForests$forest
forests <- list(forest())
for (k in c(25, 100))
  forests[[length(forests) + 1L]] <- forest(x1, basis = I(w / k))
looped <- dbarts(y ~ x1 + x3, d, forests = forests,
                 control = dbartsControl(n.chains = 1L, n.trees = 10L))
vapply(looped$data@bases[-1L], function(basis) basis[1L, 1L], 0) * c(25, 100)
#> [1] 58.93674 58.93674

# two columns, a coefficient each; I() holds arithmetic
twoColumns <- dbarts(y ~ forest(x1 + x3) + forest(x1, basis = z + I(w / 50)),
                     d, control = dbartsControl(n.chains = 1L, n.trees = 10L))
colnames(twoColumns$data@bases[[2L]])
#> [1] "z"       "I(w/50)"
```
