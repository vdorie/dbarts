# Partial Dependence Plots for BART

Run [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) at test
observations constructed so that a plot can be created displaying the
effect of a single variable (`pdbart`) or pair of variables (`pd2bart`).
Note that if \\y\\ is binary with \\P(Y=1 \| x) = F(f(x))\\, \\F\\ the
standard normal cdf, then the plots are all on the \\f\\ scale.

## Usage

``` r
pdbart(
    x.train, y.train,
    xind = NULL,
    levs = NULL, levquants = c(0.05, seq(0.1, 0.9, 0.1), 0.95),
    pl = TRUE, plquants = c(0.05, 0.95),
    ...)

# S3 method for class 'pdbart'
plot(
    x,
    xind = seq_along(x$fd),
    plquants = c(0.05, 0.95), cols = c('blue', 'black'),
    ...)

pd2bart(
    x.train, y.train,
    xind = NULL,
    levs = NULL, levquants = c(0.05, seq(0.1, 0.9, 0.1), 0.95),
    pl = TRUE, plquants = c(0.05, 0.95),
    ...)

# S3 method for class 'pd2bart'
plot(
    x,
    plquants = c(0.05, 0.95), contour.color = 'white',
    justmedian = TRUE,
    ...)
```

## Arguments

- x.train:

  Explanatory variables for training (in sample) data. Can be any valid
  input to [`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
  such as a matrix or a formula. Also accepted are fitted `bart` models
  or
  [`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  with `keepTrees` equal to `TRUE`. A sampler passed without `keepTrees`
  is used anyway, but generates fresh samples and changes its state in
  the process, with a warning; a `bart` model without a kept sampler is
  instead refit from its saved call, also with a warning (both class
  `dbartsFallbackWarning`) - and is refused outright if no call was
  saved either.

- y.train:

  Dependent variable for training (in sample) data. Can be a numeric
  vector or, when passing `x.train` as a formula, a `data.frame` or
  other object used to find variables. Not required if `x.train` is a
  fitted model or sampler.

- xind:

  Integer, character vector, or the right-hand side of a formula
  indicating which variables are to be plotted. In `pdbart`, corresponds
  to the variables (columns of `x.train`) for which a plot is to be
  constructed. In `plot.pdbart`, corresponds to the indices in list
  returned by `pdbart` for which plot is to be constructed. In
  `pd2bart`, the indices of a pair of variables (columns of `x.train`)
  to plot. If `NULL` a default of all columns is used for `pdbart` and
  the first two columns is used for `pd2bart`.

- levs:

  Gives the values of a variable at which the plot is to be constructed.
  Must be a list, where the \\i\\th component gives the values for the
  \\i\\th variable. In `pdbart`, it should have same length as `xind`.
  In `pd2bart`, it should have length 2. See also argument `levquants`.

- levquants:

  If `levs` is `NULL`, the values of each variable used in the plot are
  set to the quantiles (in `x.train`) indicated by levquants. Must be a
  vector of numeric type.

- pl:

  For `pdbart` and `pd2bart`, if `TRUE`, plot is subsequently made (by
  calling `plot.*`).

- plquants:

  In the plots, beliefs about \\f(x)\\ are indicated by plotting the
  posterior median and a lower and upper quantile. `plquants` is a
  double vector of length two giving the lower and upper quantiles.

- ...:

  Additional arguments. In `pdbart` and `pd2bart`, arguments are passed
  on to [`bart`](https://vdorie.github.io/dbarts/reference/bart.md). In
  `plot.pdbart`, they are passed on to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html). In
  `plot.pd2bart`, they are passed on to
  [`image`](https://rdrr.io/r/graphics/image.html).

- x:

  For `plot.*`, object returned from `pdbart` or `pd2bart`.

- cols:

  Vector of two colors. The first color is for the median of \\f\\,
  while the second color is for the upper and lower quantiles.

- contour.color:

  Color for contours plotted on top of the image.

- justmedian:

  A logical where if `TRUE` just one plot is created for the median of
  \\f(x)\\ draws. If `FALSE`, three plots are created one for the median
  and two additional ones for the lower and upper quantiles. In this
  case, [`mfrow`](https://rdrr.io/r/graphics/par.html) is set to
  `c(1,3)`.

## Details

We divide the predictor vector \\x\\ into a subgroup of interest,
\\x_s\\ and the complement \\x_c = x \setminus x_s\\. A prediction
\\f(x)\\ can then be written as \\f(x_s, x_c)\\. To estimate the effect
of \\x_s\\ on the prediction, Friedman suggests the partial dependence
function \$\$f_s(x_s) = \frac{1}{n}\sum\_{i=1}^n f(x_s,x\_{ic})\$\$
where \\x\_{ic}\\ is the \\i\\th observation of \\x_c\\ in the data.
Note that \\(x_s, x\_{ic})\\ will generally not be one of the observed
data points. Using BART it is straightforward to then estimate and even
obtain uncertainty bounds for \\f_s(x_s)\\. A draw of \\f^\*\_s(x_s)\\
from the induced BART posterior on \\f_s(x_s)\\ is obtained by simply
computing \\f^\*\_s(x_s)\\ as a byproduct of each MCMC draw \\f^\*\\.
The median (or average) of these MCMC draws \\f^\*\_s(x_s)\\ then yields
an estimate of \\f_s(x_s)\\, and lower and upper quantiles can be used
to obtain intervals for \\f_s(x_s)\\.

In `pdbart` \\x_s\\ consists of a single variable in \\x\\ and in
`pd2bart` it is a pair of variables.

This is a computationally intensive procedure. For example, in `pdbart`,
to compute the partial dependence plot for 5 \\x_s\\ values, we need to
compute \\f(x_s, x_c)\\ for all possible \\(x_s, x\_{ic})\\ and there
would be \\5n\\ of these where \\n\\ is the sample size. All of that
computation would be done for each kept BART draw. For this reason
running BART with `keepevery` larger than 1 (eg. 10) makes the procedure
much faster.

## Value

The plot methods produce the plots and don't return anything.

`pdbart` and `pd2bart` return lists with components given below. The
list returned by `pdbart` is assigned class `pdbart` and the list
returned by `pd2bart` is assigned class `pd2bart`. The class carries a
`plot` method only; it is not a fit, so `predict`, `extract`, `fitted`,
and `residuals` are not defined for it (`fitted` and `residuals` fall
through to stats' defaults and return `NULL`).

- fd:

  A matrix whose \\(i, j)\\ value is the \\i\\th draw of \\f_s(x_s)\\
  for the \\j\\th value of \\x_s\\. “fd” is for “function draws”.

  For `pdbart` `fd` is actually a list whose \\k\\th component is the
  matrix described above corresponding to the \\k\\th variable chosen by
  argument `xind`. The number of columns in each matrix will equal the
  number of values given in the corresponding component of argument
  `levs` (or number of values in `levquants`).

  For `pd2bart`, `fd` is a single matrix. The columns correspond to all
  possible pairs of values for the pair of variables indicated by
  `xind`. That is, all possible \\(x_i, x_j)\\ where \\x_i\\ is a value
  in the levs component corresponding to the first \\x\\ and \\x_j\\ is
  a value in the levs components corresponding to the second one. The
  first \\x\\ changes first.

- levs:

  The list of levels used, each component corresponding to a variable.
  If argument `levs` was supplied it is unchanged. Otherwise, the levels
  in `levs` are as constructed using argument `levquants`.

- xlbs:

  A vector of character strings which are the plotting labels used for
  the variables.

The remaining components are passed on from the BART run used to create
the partial dependence plot, under
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s own names
(see its ‘Value’ section) with one rename and several omissions - the
fit's matched call is under `bartcall` rather than `call`, and not every
`bart` component is carried over. For a continuous response the
remaining components are `bartcall`, `yhat.train`, `first.sigma`,
`sigma`, `yhat.train.mean`, `sigest`, `y`, and `fit`; `call`, `family`,
`resid.dist`, `yhat.test`, `yhat.test.mean`, `s.train`, `s.test`,
`varcount`, and `n.chains` are not present. For a binary response only
`bartcall`, `yhat.train`, `y`, and `fit` are present; `first.sigma`,
`sigma`, `yhat.train.mean`, `sigest`, and `binaryOffset` are absent as
well, on top of the continuous-case omissions. `pd2bart` follows the
same rules. The function
[`plot.bart`](https://vdorie.github.io/dbarts/reference/bart.md) can be
applied to the object returned by `pdbart` or `pd2bart` to examine the
BART run.

## References

Chipman, H., George, E., and McCulloch, R. (2010) BART: Bayesian
additive regression trees. *The Annals of Applied Statistics*, **4**(1),
266–298. [doi:10.1214/09-AOAS285](https://doi.org/10.1214/09-AOAS285) .

Chipman, H., George, E., and McCulloch R. (2006) Bayesian Ensemble
Learning. Advances in Neural Information Processing Systems 19,
Scholkopf, Platt and Hoffman, Eds., MIT Press, Cambridge, MA, 265-272.
<https://www.rob-mcculloch.org/>

Friedman, J.H. (2001) Greedy function approximation: A gradient boosting
machine. *The Annals of Statistics*, **29**, 1189–1232.

## Author

Hugh Chipman: <hugh.chipman@gmail.com>.  
Robert McCulloch: <robert.mcculloch1@gmail.com>.

## Examples

``` r
## simulate data
f <- function(x)
    return(0.5 * x[,1] + 2 * x[,2] * x[,3])

sigma <- 0.2
n     <- 100

set.seed(27)
x <- matrix(2 * runif(n * 3) - 1, ncol = 3)
colnames(x) <- c('rob', 'hugh', 'ed')

Ey <- f(x)
y  <- rnorm(n, Ey, sigma)

# \donttest{
## pdbart: one dimensional partial dependence plot
set.seed(99)
pdb1 <- pdbart(
    x, y, xind = c(1, 2),
    levs = list(seq(-1, 1, 0.2), seq(-1, 1, 0.2)),
    pl = FALSE, keepevery = 10, ntree = 100, verbose = FALSE
)

## pd2bart: two dimensional partial dependence plot
set.seed(99)
pdb2 <- pd2bart(
    x, y, xind = c(2, 3),
    levquants = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
    pl = FALSE, ntree = 100, keepevery = 10, verbose = FALSE)

## the plot methods draw into the current device, so set up the layout and
## put the caller's graphical parameters back afterwards. The first two plot
## regions are for pdbart, the third for pd2bart.
oldpar <- par(no.readonly = TRUE)
par(mfrow = c(1, 3))
plot(pdb1, ylim = c(-0.6, 0.6))
plot(pdb2)

par(oldpar)

## compare BART fit to linear model and truth = Ey
lmFit  <- lm(y ~ ., data.frame(x, y))
fitmat <- cbind(y, Ey, lmFit$fitted, pdb1$yhat.train.mean)
colnames(fitmat) <- c('y', 'Ey', 'lm', 'bart')
print(cor(fitmat))
#>              y        Ey        lm      bart
#> y    1.0000000 0.9603886 0.4052732 0.9900234
#> Ey   0.9603886 1.0000000 0.4457354 0.9813192
#> lm   0.4052732 0.4457354 1.0000000 0.4375004
#> bart 0.9900234 0.9813192 0.4375004 1.0000000
# }

# \donttest{
## example showing the use of a pre-fitted model
df <- data.frame(y, x)
set.seed(99)
bartFit <- bart(
    y ~ rob + hugh + ed, df,
    keepevery = 10, ntree = 100, keeptrees = TRUE, verbose = FALSE)
pdb3 <- pdbart(bartFit, xind = rob + ed, pl = FALSE)
# }
```
