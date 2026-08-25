# Discrete Bayesian Additive Regression Trees Sampler Data

Convenience function to create a data object for use with a
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) sampler.

## Usage

``` r
dbartsData(
    formula, data, test, subset, weights,
    offset, offset.test = offset,
    factors = c("categorical", "indicators"),
    missing = c("incorporate", "error"),
    bases = NULL,
    counts = NULL)
```

## Arguments

- formula, data, test, subset, weights, offset, offset.test, factors,
  missing:

  As in [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md).
  Retains backwards compatibility with
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md), so that
  `formula`/`data` can be a
  [`formula`](https://rdrr.io/r/stats/formula.html)/[`data.frame`](https://rdrr.io/r/base/data.frame.html)
  pair, or a pair of `x.train`/`y.train` matrices/vector.
  `offset`/`offset.test` take the matrix interface's shape from whether
  `counts` is also supplied: with `counts`, an `n` x `K` (`nTest` x `K`)
  numeric matrix installs the per-category shift described under
  `counts` below, and a flat vector is refused; without `counts`, only a
  flat vector (or single number) is accepted, and a matrix is refused.

- bases:

  Optional list of the per-forest amplitude bases a multi-forest model
  combines its forests through, one element per forest and `NULL` where
  a forest carries the implicit intercept its single amplitude scales.
  Each non-null element is a numeric matrix with one row per observation
  and one amplitude per column, restricted by `subset` exactly as
  `weights` is. This is the data half of a multi-forest specification,
  and is normally written by
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
  `forests = ` argument rather than here; supplying it directly selects
  the multi-forest model with the engine's default per-forest
  configuration, and
  [`forest`](https://vdorie.github.io/dbarts/reference/forest.md)
  declarations may then configure the forests one at a time.

- counts:

  Optional `n` x `K` matrix of non-negative integer counts selecting the
  multinomial (softmax) model, column `k` holding category `k`'s
  successes and every row carrying at least one trial. It is restricted
  by `subset` exactly as `weights` is, and it IS the response: `y` is
  derived as its row sums - the trials \\n_i\\ - so no separate response
  need be supplied, and any that is supplied is replaced. At least two
  categories are required. Column names, when present, label the
  categories.
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)`(family = "multinomial")`
  writes this argument from a count-matrix or factor response rather
  than requiring it here. Paired with `counts`, a matrix
  `offset`/`offset.test` is the per-category shift added to the raw
  per-category values BEFORE the softmax - which is the only place a
  shift means anything under one, since the blend is invariant to a
  common per-observation shift and a flat shift would be added after it.
  Only the row-centred part is identified: adding a constant to a whole
  row leaves every reported probability unchanged.

## Details

### Data frame ingestion

Data frame columns - through the formula interface or a frame passed as
`x.train` - map to predictors as follows:

- numeric and logical columns enter as ordinal predictors;

- unordered factors each become a single categorical predictor under the
  default `factors = "categorical"`, with splits sending a subset of
  levels down each branch; `factors = "indicators"` instead expands each
  factor into binary indicator columns, as previous versions always did;

- ordered factors enter as ordinal predictors, coded numerically by
  level order;

- matrix columns (e.g. from [`poly`](https://rdrr.io/r/stats/poly.html))
  splice in as ordinal predictors, one per column;

- [`Matrix::sparseVector`](https://rdrr.io/pkg/Matrix/man/sparseVector.html)
  and `dgCMatrix` columns (assigned into the frame; they do not survive
  `data.frame(...)` or [`I()`](https://rdrr.io/r/base/AsIs.html)) enter
  as sparse ordinal predictors. Only the x/y interface accepts them for
  training: a bare S4 column does not survive
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html), so the
  formula path refuses one explicitly, the same guarded message
  [`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)
  columns get (see below), rather than raw-erroring inside
  `model.frame`. `test`/`x.test` accepts them under either interface,
  since a test set is never run through `model.frame`;

- [`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)
  columns enter as a single sparse categorical predictor, coded and
  binned identically to a dense factor of the same values. Only the x/y
  interface accepts them (a bare S4 column cannot survive
  [`model.frame`](https://rdrr.io/r/stats/model.frame.html), so the
  formula path refuses one explicitly).

Factor level tables are retained so that test data are coded
identically, and a categorical predictor's category count is the
*declared* level table rather than the levels the training rows happen
to take. A level declared but unobserved in training therefore keeps its
own bin, and test or replacement data carrying it are accepted at every
entrance (creation, `setTestPredictor`, `predict`, `setPredictor`); this
holds equally for dense factors, `sparseFactor` columns, and dense
factors inside a mixed dense/sparse frame. Drop unused levels before
fitting ([`droplevels`](https://rdrr.io/r/base/droplevels.html)) if you
would rather they not be modeled. `test` accepts the same column types
as `x.train`, including `sparseFactor` and sparse ordinal columns. A
`sparseFactor` or sparse ordinal test column recodes over the training
level table and stays resident - through creation and
`setTestPredictor` - rather than densifying at ingestion, the same
storage tier rule as training. `predict` and `getTrees(newdata = )` code
a test set the same way and then materialize it to a dense matrix at the
call boundary, since those entry points evaluate already-built trees on
resolved values.

Data frame input is stored columnar: no \\n \times p\\ double matrix is
retained. Code that previously reached into the data object's `x` slot
expecting a plain matrix should call
[`extract`](https://vdorie.github.io/dbarts/reference/bart.md)`(sampler, "predictors")` -
or `as.matrix` on the slot - to obtain the numeric predictor matrix,
factor columns as their integer codes.

## Value

An object of class `dbartsData`.

## See also

[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`extract.dbartsSampler`](https://vdorie.github.io/dbarts/reference/extract.dbartsSampler.md),
[`sparseFactor`](https://vdorie.github.io/dbarts/reference/sparseFactor.md)

## Examples

``` r
set.seed(3)
n <- 60L
df <- data.frame(
  x = runif(n),
  g = factor(sample(c("a", "b", "c"), n, replace = TRUE))
)
df$y <- df$x + 0.5 * (df$g == "b") + rnorm(n, 0, 0.1)

## the formula interface; 'g' enters as one categorical predictor
data <- dbartsData(y ~ x + g, df)
data
#> An object of class "dbartsData"
#> Slot "y":
#>  [1]  0.177577428  0.760531131  0.993684793  0.265626355  0.697894896
#>  [6]  0.519816750  0.009659326  0.244830643  1.014011813  1.005901046
#> [11]  0.335204293  0.506964701  0.651273625  1.127721920  1.276516351
#> [16]  0.841611818 -0.011854548  1.330332559  0.905008414  0.857943974
#> [21]  0.244923711 -0.034990768  0.499862880  0.052575161  0.621328776
#> [26]  0.745552247  0.703889559  1.424777615  1.032684790  1.387674186
#> [31]  0.320629695  0.481533368  0.168516403  0.924658735  0.806339550
#> [36]  0.152153052  0.883719720  0.120276483  1.272290346  0.589972076
#> [41]  0.344054654  0.874243472  0.697341157  0.628454551  0.279732169
#> [46]  0.232456352  0.026730705  0.133342322  0.441749469  1.251978519
#> [51]  0.253527651  0.275760383  0.727522408  1.454490179  1.054661032
#> [56]  0.964278766  0.405456058  0.254249106  0.132526675  0.272371030
#> 
#> Slot "x":
#> $dense
#> $dense[[1]]
#>  [1] 0.16804153 0.80751640 0.38494235 0.32773432 0.60210067 0.60439405
#>  [7] 0.12463344 0.29460092 0.57760992 0.63097927 0.51201590 0.50502391
#> [13] 0.53403535 0.55724944 0.86791949 0.82970869 0.11144915 0.70368836
#> [19] 0.89748826 0.27973255 0.22820188 0.01532989 0.12898156 0.09338193
#> [25] 0.23688501 0.79114741 0.59973157 0.91014771 0.56042455 0.75570477
#> [31] 0.37917189 0.37328098 0.17029064 0.45330732 0.25841396 0.33626595
#> [37] 0.88958303 0.20194630 0.57918604 0.20763203 0.28146879 0.78628120
#> [43] 0.17301935 0.57074752 0.41928296 0.26762217 0.04780944 0.10349305
#> [49] 0.31403146 0.80064106 0.22932470 0.21299844 0.87710091 0.99322196
#> [55] 0.84424702 0.91043655 0.47126973 0.22441841 0.12781466 0.27968351
#> 
#> $dense[[2]]
#>  [1] a c b a a a a a b b a c c b b a c b a b a a b a b c c b b b c c a b b a c c
#> [39] b b a a b c c c c a a b c a c b a c c a c a
#> Levels: a b c
#> 
#> 
#> $sparse
#> NULL
#> 
#> $map
#> [1] 1 2
#> 
#> $numObservations
#> [1] 60
#> 
#> $columnNames
#> [1] "x" "g"
#> 
#> attr(,"class")
#> [1] "dbartsMixedMatrix"
#> attr(,"term.labels")
#> [1] "x" "g"
#> attr(,"varTypes")
#> [1] 0 1
#> attr(,"factor.levels")
#> attr(,"factor.levels")[[1]]
#> NULL
#> 
#> attr(,"factor.levels")[[2]]
#> [1] "a" "b" "c"
#> 
#> 
#> Slot "varTypes":
#> [1] 0 1
#> 
#> Slot "x.test":
#> NULL
#> 
#> Slot "weights":
#> NULL
#> 
#> Slot "weights.test":
#> NULL
#> 
#> Slot "offset":
#> NULL
#> 
#> Slot "offset.test":
#> NULL
#> 
#> Slot "n.cuts":
#> [1] NA NA
#> 
#> Slot "sigma":
#> [1] NA
#> 
#> Slot "missing":
#> [1] "incorporate"
#> 
#> Slot "response.type":
#> [1] "numeric"
#> 
#> Slot "response.n.levels":
#> [1] NA
#> 
#> Slot "response.levels":
#> NULL
#> 
#> Slot "bases":
#> NULL
#> 
#> Slot "counts":
#> NULL
#> 
#> Slot "offset.category":
#> NULL
#> 
#> Slot "offset.category.test":
#> NULL
#> 
#> Slot "testUsesRegularOffset":
#> [1] NA
#> 

## the same data through the backwards-compatible x/y interface, with a
## test set that is coded against the training level table
train <- df[1:40, ]
test  <- df[41:60, c("x", "g")]
data <- dbartsData(train[, c("x", "g")], train$y, test)

## build a sampler on it, rather than letting dbarts() ingest again
sampler <- dbarts(data, control = dbartsControl(
  n.chains = 1L, n.threads = 1L, n.samples = 10L, n.burn = 10L,
  n.trees = 25L, updateState = FALSE))
samples <- sampler$run()
str(samples$test)
#>  num [1:20, 1:10] 0.347 0.933 0.469 0.612 0.39 ...
```
