# R-hat of One Quantity's Draws

The rank-normalized split R-hat (Vehtari, Gelman, Simpson, Carpenter and
Buerkner, 2021) of the draws of one quantity from several chains: the
larger of the bulk and the folded (tail) values. It is the statistic
[`summary.bart`](https://vdorie.github.io/dbarts/reference/summary.bart.md)
reports, and from four iterations on it returns the value
`posterior::rhat` does for the inputs that function takes (see Details
for the inputs where the two differ).

## Usage

``` r
rhat(x)
```

## Arguments

- x:

  A numeric (or logical) matrix of the draws of one quantity, with a row
  for each iteration and a column for each chain. A vector is taken as
  one chain.

## Details

Each chain is cut in two, the middle draw of an odd length dropped, and
the draws are rank-normalized as one pool; the bulk value compares the
half-chains' means and variances, and the tail value does the same on
the draws folded about their median. R-hat compares each chain's spread
with the pooled spread, so a value above 1 says the chains disagree.

The result is `NA`, as for `posterior::rhat`, when `x` has an `NA`,
`NaN` or constant draws, or has a single iteration. A non-finite draw
other than `NA` and `NaN` is ranked, as posterior ranks it. An argument
that is not numeric or logical, is an array of more than two dimensions,
or has no iteration or no chain is an error. posterior refuses a data
frame and a character vector too, but takes a complex vector (and
returns an infinite value with a warning for a matrix with no chain);
dbarts refuses those, as they are not draws.

With two or three iterations and several chains, posterior (1.7.0) loses
the dimension of its one-draw half-chains and treats each iteration as a
chain; dbarts returns `NA`, as a half-chain of one draw has no variance.
At four or more iterations the two agree to rounding.

R-hat should be read on the estimand reported: at the defaults BART's
chains routinely disagree on `sigma`, `k` and single predictions, while
averaged estimands, such as a subgroup's mean prediction or a contrast
between subgroups, mix far better. See
[`summary.bart`](https://vdorie.github.io/dbarts/reference/summary.bart.md).

## Calling `rhat`

Call it as `dbarts::rhat`. Other packages export a function named
`rhat`: posterior and bayesplot (and brms re-exports posterior's), and
rstan exports `Rhat`. Whichever is attached last masks the others, and
bayesplot's takes a fitted Stan model, not a matrix. The stan4bart
package attaches dbarts, so it brings this name in.

## Value

A single number, `NA_real_` where the draws carry no information.

## References

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Buerkner, P.-C.
(2021) Rank-normalization, folding, and localization: an improved
\\\hat{R}\\ for assessing convergence of MCMC. *Bayesian Analysis*
16(2), 667-718.

## See also

[`summary.bart`](https://vdorie.github.io/dbarts/reference/summary.bart.md),
[`extract`](https://vdorie.github.io/dbarts/reference/bartBT.md)

## Examples

``` r
set.seed(1)
draws <- matrix(rnorm(400), 100L, 4L)
dbarts::rhat(draws)
#> [1] 0.9948036

# four chains that disagree
dbarts::rhat(draws + rep(0:3, each = 100L))
#> [1] 1.528219
```
