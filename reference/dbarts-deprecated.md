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
error.

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

- `resid.dist`, `dispersion`, `breaks`, `max.rows`:

  On [`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md), and
  (the first two)
  [`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md):
  family-only settings, now arguments of the
  [`dbartsFamilies`](https://vdorie.github.io/dbarts/reference/dbartsFamilies.md)
  constructor that owns each (`family = student(df)`,
  `family = nbinom(dispersion)`, `family = hazard(breaks)`,
  `family = hazard(max.rows)`).

- `dart`, `levelGibbs`:

  On [`bart`](https://vdorie.github.io/dbarts/reference/bart.md) (`dart`
  also on
  [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md)):
  tree-prior-only settings, now arguments of
  [`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)'s
  tree priors (`tree.prior = dart()`,
  `tree.prior = cgm(levelGibbs = )`).

- `control` on
  [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md):

  Refused outright, naming the flat arguments its settings became
  (`n.cuts`, `useQuantiles`, `n.thin`, `storage`): `xbart` builds its
  own control, so there is nothing left for one to set.

- A three-element `n.burn` on
  [`xbart`](https://vdorie.github.io/dbarts/reference/xbart.md):

  Refused outright, naming `n.burn = c(fresh, warm)`: 0.9-x's third
  element was a per-replication burn-in, and a chain is never carried
  between replications now.

- `family = "twopart"`:

  Refused outright on every entry point that takes `family`, naming
  `family = "hurdle.lognormal"`.

- `dbartsSampler`'s `$startThreads`/`$stopThreads`:

  No-op methods
  ([`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)):
  the hierarchical thread manager they drove is gone, threads now owned
  per run, but a 0.9-x Gibbs loop that brackets its sweeps with them
  still runs.

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

Every entry above, and `bart2` and the `rbart_vi` stubs, is one row of
`dbarts:::dbartsTombstones`, the registry this page documents by hand.

## Arguments

- object:

  An `rbart` fit object saved by dbarts 0.9-x.

- ...:

  Ignored; the stubs do not run. `bart2`'s own arguments are
  [`bart`](https://vdorie.github.io/dbarts/reference/bart.md)'s, copied
  rather than restated - `formals(bart2)` is `formals(bart)` - so a
  caller reading the defaults programmatically still finds them.

## Value

`bart2` returns whatever
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md) returns. The
`rbart_vi` stubs do not return.

## See also

[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md)
