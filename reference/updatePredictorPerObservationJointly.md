# Jointly Update a Shared Predictor per Observation Across Samplers

Applies a “partial” (per-observation) update of a single shared column
to several
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
objects at once, installing an observation's new value only when it
keeps every leaf node non-empty in every tree of every sampler.

## Usage

``` r
updatePredictorPerObservationJointly(samplers, x, column, updateState = NULL)
```

## Arguments

- samplers:

  A
  [`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md),
  or a list of `dbartsSampler` objects, that share the same
  (index-aligned) observations. The samplers are treated as peers; none
  is privileged.

- x:

  The new values for the shared column, one per observation. For a
  numeric column, numbers; a factor, or text that is not numerals, is
  refused for it. For a column the samplers hold as a factor or an
  ordered factor, its labels: a factor, a character vector or a
  `sparseFactor`, matched to the column's levels by name, as
  [`setPredictor`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  takes them for a named column. The order of the levels of the factor
  given does not matter; a label the column does not have is refused by
  name, and so is a missing value when the column holds none. Such a
  column also takes numbers, read as the codes each sampler holds it in:
  whole numbers from 0 in the order of the column's levels, so the first
  level is 0. `as.integer(f)` of a factor counts from 1 and is not these
  codes, so give the factor itself. A number that is not a code is
  refused, and a missing number is installed as a missing value, where a
  missing label is refused in a column that has none. A number is never
  matched to a label, even where the labels are numerals, so give labels
  as a factor or as character. One value is one level in every sampler
  only if they hold the column with the same levels in the same order,
  so samplers that hold it as a factor with different levels are
  refused, for labels and numbers alike; update them in separate calls.
  Labels are also refused when one sampler holds the column as a number,
  where numbers are each sampler's own.

- column:

  A single integer or character string identifying the shared column; it
  cannot be missing. The column is matched *by name* across all
  samplers, so the shared variable may sit at different column positions
  in each. When given as an integer, it indexes the columns of the first
  sampler and the corresponding name is used to locate the column in the
  rest.

- updateState:

  A logical determining if the local cache of each sampler's state
  should be updated after the update completes. If `NULL` (the default),
  each sampler's own
  [`control`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  default is used; `TRUE` or `FALSE` overrides it and anything else is
  an error. An explicit `NA` reads as `NULL` for one release, with a
  once-per-session warning.

## Details

This is the multi-sampler analog of the `forceUpdate = "partial"` mode
of
[`setPredictor`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md).
A single sequential sweep installs each observation in every sampler at
once, and only if its new value keeps every leaf non-empty in every tree
of every *forest* of every chain of every sampler; otherwise that
observation is rolled back to its previous value in all of them. In a
sampler with a
[`monotone`](https://vdorie.github.io/dbarts/reference/monotone.md)
constraint a value is also declined, in every sampler, when installing
it would leave a tree's leaf values out of order; that page says which
values do this and why proposing them again may never succeed. Any of
the samplers may itself carry more than one forest - a Bayesian causal
forest, a multinomial sampler, or a heteroscedastic (`variance`)
sampler - and each of its forests is checked on the same terms as an
ordinary single-forest sampler's mean forest (see ‘Multi-forest and
heteroscedastic predictor mutation’ in
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)).
Like the single-sampler partial mode, the update never changes tree
structure – it only re-routes observations.

This supports embedding several BART components that condition on a
common latent predictor (for example, a shared ability or trait in an
item-response model) within a larger Gibbs/Metropolis sampler.

A whole-matrix `setPredictor` call on one of the samplers, made earlier
in the same script, keeps that sampler's column names (see ‘Multi-forest
and heteroscedastic predictor mutation’,
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)),
so this function's by-name `column` matching continues to resolve
against it afterward.

## Value

A logical vector of length equal to the number of observations, `TRUE`
where the observation's new value was installed and `FALSE` where it was
rolled back to keep every tree valid. The same accept/reject decision is
applied to every sampler, so an installed observation takes the new
value in all of them; a rejected observation reverts to each sampler's
own previous value. When the samplers agreed on the shared column before
the call, they therefore continue to agree afterward.

## See also

[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)

## Examples

``` r
## two BART components conditioning on a common latent predictor 'theta',
## as an item-response model would: both must accept an observation's new
## value, or neither does
set.seed(5)
n <- 80L
theta <- rnorm(n)
w <- runif(n)

df1 <- data.frame(theta = theta, w = w)
df1$y <- theta + rnorm(n, 0, 0.2)
df2 <- data.frame(v = runif(n), theta = theta)
df2$y <- 0.5 * theta + rnorm(n, 0, 0.2)

control <- dbartsControl(n.chains = 1L, n.threads = 1L, n.burn = 0L,
                         n.samples = 1L, n.trees = 25L, updateState = FALSE)
s1 <- dbarts(y ~ theta + w, df1, control = control)
s2 <- dbarts(y ~ v + theta, df2, control = control)
invisible(s1$run())
invisible(s2$run())

## 'theta' sits in a different column of each sampler; it is matched by name
proposal <- theta + rnorm(n, 0, 0.1)
installed <- updatePredictorPerObservationJointly(
  list(s1, s2), proposal, column = "theta")
table(installed)
#> installed
#> TRUE 
#>   80 
```
