# Block-Additive Constraints for BART

Build a block-additive constraint confining each whole tree in a BART
fit to one declared group of predictors, so the ensemble is exactly a
sum of per-group functions \\f = \sum_G f_G\\ (a functional-ANOVA /
grouped-GAMI decomposition). Pass the result as the `blocks` argument of
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) or
[`bart2`](https://vdorie.github.io/dbarts/reference/bart.md). The
constraint is applied per forest, so in a multi-forest model each forest
can carry its own partition.

## Usage

``` r
blocks(groups, trees.per.group = NULL)
```

## Arguments

- groups:

  A list partitioning the forest's predictors into disjoint blocks. Each
  entry is a character vector of model-matrix column names (a bare
  factor term name expands to its indicator columns) or a numeric vector
  of column indices. The partition must be *total* and *disjoint*: every
  predictor is named in exactly one group. A single vector (rather than
  a list) is taken as one group.

- trees.per.group:

  Optional integer vector, one entry per group, giving how many of the
  `n.trees` trees each group receives; the entries must be positive and
  sum to `n.trees`. `NULL` (the default) distributes the trees as evenly
  as possible, with the first `n.trees %% length(groups)` groups
  receiving one extra tree.

## Details

Each tree is assigned to exactly one group and split-restricted to that
group's predictors, so no path in the tree - and hence no term in the
ensemble - mixes predictors from different groups. The assignment is
deterministic and contiguous (the first `trees.per.group[1]` trees to
group 1, the next `trees.per.group[2]` to group 2, and so on) and
consumes no random numbers, so a fit is exactly reproducible.

Unlike
[`interactions`](https://vdorie.github.io/dbarts/reference/interactions.md)`(groups = )`,
which is a per-path allow-list that leaves un-named predictors
unconstrained, `blocks` lowers to a static per-tree column mask: a
predictor named in no block would be masked out of every tree and never
enter the fit. The partition is therefore required to be total, and
every value is validated against the model matrix at fit time - an
unpartitioned predictor, a predictor named in more than one group, an
unrecognized name, or a `trees.per.group` of the wrong length, with a
non-positive entry, or not summing to `n.trees` is an error. For a
forest that is already restricted to a subset of the columns, the
partition is over that forest's available columns.

`blocks` and
[`interactions`](https://vdorie.github.io/dbarts/reference/interactions.md)
act on independent seams (which columns a tree may use versus which may
co-occur on a path) and may be combined on one forest. Block-additivity
is also orthogonal to monotone constraints (see
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)'s
`monotone`).

## Value

A `dbartsBlocks` specification object, resolved when a sampler is built.

## See also

[`interactions`](https://vdorie.github.io/dbarts/reference/interactions.md),
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
[`bart2`](https://vdorie.github.io/dbarts/reference/bart.md)

## Examples

``` r
set.seed(0)
n <- 100L
x <- matrix(runif(n * 3), n, 3, dimnames = list(NULL, c("x1", "x2", "x3")))
y <- 2 * x[, 1] + ifelse(x[, 2] > 0.5, 1, -1) + rnorm(n, 0, 0.2)
df <- data.frame(y, x)

## block-additive fit: f(x) = f1(x1) + f2(x2, x3), no tree mixing the blocks
fit <- bart2(y ~ x1 + x2 + x3, df,
             blocks = blocks(groups = list("x1", c("x2", "x3"))),
             n.trees = 20L, n.samples = 20L, n.burn = 20L,
             n.chains = 1L, verbose = FALSE)

## split the 20 trees explicitly: 5 to the first block, 15 to the second
fit.split <- bart2(y ~ x1 + x2 + x3, df,
                   blocks = blocks(groups = list("x1", c("x2", "x3")),
                                   trees.per.group = c(5L, 15L)),
                   n.trees = 20L, n.samples = 20L, n.burn = 20L,
                   n.chains = 1L, verbose = FALSE)
```
