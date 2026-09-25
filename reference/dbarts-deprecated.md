# Deprecated and Removed Functions in dbarts

Names kept reachable for one release past their removal or rename, each
pointing at its successor. They are removed in dbarts 1.1-0.

`bart2` is the former name of the modern front door, which is now
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md); a `bart2`
call is forwarded unchanged, with a warning the first time in a session.

`rbart_vi` and its methods fit a grouped random-intercept BART model.
That model now lives in the stan4bart package (`stan4bart::stan4bart`),
whose prior on the group spread differs from the one dbarts drew from,
so a refit moves the results rather than reproducing them. These stubs
error, except `print`, which prints the stored call, as 0.9-x did, and a
note naming stan4bart.

## Usage

``` r
rbart_vi(...)

# S3 method for class 'rbart'
predict(object, ...)
# S3 method for class 'rbart'
extract(object, ...)
# S3 method for class 'rbart'
fitted(object, ...)
# S3 method for class 'rbart'
residuals(object, ...)
# S3 method for class 'rbart'
plot(x, ...)
# S3 method for class 'rbart'
print(x, ...)
```

## Details

The rest of the registry is names still reachable on entry points that
keep their own page, rather than functions of their own: each is
accepted for one release, mapped onto its successor with a
once-per-session warning (or refused by name, where noted), and removed
in dbarts 1.1-0.

- `rngSeed`:

  On [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) and
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md):
  the engine seed, now `seed`.

- `sigma`:

  On [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) and
  [`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md):
  the creation-time residual estimate, now `sigest`.

- `power`, `base`:

  On [`bart`](https://vdorie.github.io/dbarts/reference/bart.md):
  tree-prior-only settings, now arguments of
  [`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)'s
  `cgm()` or `dart()` (`tree.prior = cgm(power, base)` or
  `tree.prior = dart(power, base)`).

- `split.probs`:

  On [`bart`](https://vdorie.github.io/dbarts/reference/bart.md): a
  tree-prior-only setting, now an argument of
  [`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)'s
  `cgm()` only (`tree.prior = cgm(split.probs = )`); a DART tree prior
  draws its own split probabilities, so combining it with
  `tree.prior = dart()` is refused.

- `resid.prior`, `sigdf`, `sigquant`:

  The residual scale's prior is a setting of the families that draw one,
  so it rides the family object
  (`family = gaussian(sigma = chisq(df, quant))`,
  `family = gaussian(sigma = fixed(value))`); likewise
  `student(df, sigma = )`, `aft(sigma = )` and
  `hurdle.lognormal(sigma = )`. `sigdf` and `sigquant` are retired on
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md);
  `resid.prior` is retired the same way on
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md),
  `dbartsSpec` and
  [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md), so the
  prior has one home on every door. A call that writes it both ways - a
  retired spelling and a family whose call named `sigma` - is refused
  where the two disagree. `sigest`, the residual-scale estimate supplied
  at creation, is unaffected.

- `proposal.probs`:

  On [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) and
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md): the
  tree-move mixture, now a
  [`dbartsControl`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  slot and formal (`control = dbartsControl(proposal.probs = )`).

- A three-element `n.burn` on
  [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md):

  Refused outright, naming `n.burn = c(fresh, warm)`: 0.9-x's third
  element was a per-replication burn-in, and a chain is never carried
  between replications now.

- `dbartsSampler`'s `$startThreads`/`$stopThreads`:

  No-op methods
  ([`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)):
  the hierarchical thread manager they drove is gone, threads now owned
  per run, but a 0.9-x Gibbs loop that brackets its sweeps with them
  still runs.

- A thread count on `dbartsSampler`'s `$run`:

  Ignored after a once-per-session warning, whether passed as
  `n.threads`, as 0.9-x's formal `numThreads`, or positionally after
  `updateState`: a run uses the sampler's own count, which `$setControl`
  changes, and its draws do not depend on it.

- A BayesTree-spelled
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) call:

  Forwarded to
  [`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md) after
  a once-per-session warning, when the call carries a name from
  `bartBT`'s own formals that `bart` does not take.

- A fourth or later positional
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) argument:

  Refused, naming both doors: 0.9-x's `bart` read a fourth positional
  argument as `sigest`, where `bart`'s own fourth formal is `subset`.

- The front-door startup message:

  Attaching the package prints one `packageStartupMessage` naming both
  doors, once per session (`suppressPackageStartupMessages` silences
  it); no successor, since nothing replaces it.

- The front-door defaults message:

  A [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) call
  that names no `bartBT` argument and whose first argument is not a
  formula - a 0.9-x call such as `bart(x.train, y.train, x.test)` -
  binds the same way under both doors, so nothing forwards it and it
  runs under `bart`'s own defaults. The first such call in a session
  prints a message of class `dbartsFrontDoorMessage` naming the change
  and [`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md);
  `suppressMessages` silences it, and a call made from package code
  prints nothing.

`control` on
[`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md) is no
longer a tombstone: both `bart` and `xbart` now take a
`control = `[`dbartsControl()`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
formal, reaching every setting neither door spells flatly, so the
earlier retirement (which refused the argument outright) is reversed.

Every entry above, and `bart2` and the `rbart_vi` stubs, is one row of
`dbarts:::dbartsTombstones`, the registry this page documents by hand.

## Arguments

- object, x:

  An `rbart` fit object saved by dbarts 0.9-x.

- ...:

  Ignored; the stubs do not run. `bart2`'s own arguments are
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s, copied
  rather than restated - `formals(bart2)` is `formals(bart)` - so a
  caller reading the defaults programmatically still finds them.

## Value

`bart2` returns whatever
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md) returns. The
`rbart_vi` stubs do not return, except `print`, which returns `x`
invisibly.

## See also

[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md)
