# Augmentation Draws For A Composed Sampler

The per-observation data augmentation a non-gaussian
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
runs inside each sweep, callable from R against R's own random number
stream. `dbartsDrawLatents` draws the augmentation variable at a
supplied fit - the quantity
[`getLatents`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
reports for that family - and `dbartsWorkingResponse` turns a drawn
latent into the response a host regresses on. Neither touches a sampler:
a host writing its own Gibbs step around BART, or around anything else,
need not reimplement the family's augmentation.

## Usage

``` r
dbartsDrawLatents(family, fit, y, weights = NULL, offset = NULL, sigma = 1,
                  dispersion = NULL, cutpoints = NULL, df = NULL)

dbartsWorkingResponse(family, latent, y, weights = NULL, offset = NULL,
                      dispersion = NULL)
```

## Arguments

- family:

  One of `"probit"`, `"logistic"`, `"ordinal"`, `"aft"`, `"nbinom"` or
  `"student"`, matched by
  [`match.arg`](https://rdrr.io/r/base/match.arg.html); there is no
  default. `"student"` is the Student-t residual distribution
  (`resid.dist = student()` in
  [`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md)), a
  gaussian response under a scale mixture.

- fit:

  A numeric vector of each observation's location *without* the offset -
  \\f(x_i)\\, what
  [`getFitsWithoutOffset`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  reports, or however the host computes it. Its length sets \\n\\.

- latent:

  A numeric vector of length \\n\\ of drawn latents, as
  `dbartsDrawLatents` returns; positive for `"logistic"` and `"nbinom"`,
  whose latents are precisions.

- y:

  The response, of length \\n\\, constrained exactly as the sampler
  constrains it: 0 or 1 for `"probit"` and `"logistic"`, an integer
  category index in \\\[1, K\]\\ for `"ordinal"`, a non-negative integer
  count for `"nbinom"`. For `"aft"` it is the log observed (censoring)
  time, the bound the drawn log survival time is truncated below at;
  `"student"` constrains nothing.

- weights:

  `"logistic"` only, refused by name elsewhere: observation counts,
  positive whole numbers, each entering as the sum of that many
  independent \\PG(1, \psi)\\ draws.

- offset:

  A numeric vector of length \\n\\ added to `fit` to form the linear
  predictor, or `NULL` for none.

- sigma:

  `"aft"` and `"student"` only: the residual standard deviation, in
  response units.

- dispersion:

  `"nbinom"` only, and required there: the current dispersion \\r\\, a
  positive whole number, as
  [`getDispersion`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
  reports it.

- cutpoints:

  `"ordinal"` only, and required there: the \\K - 1\\ strictly
  increasing category boundaries in force.

- df:

  `"student"` only, and required there: the degrees of freedom \\\nu\\.

## Details

**The offset convention.** `fit` is the location *without* the offset,
the engine's own convention at every augmentation site, and the linear
predictor is formed internally as `fit + offset`. The two arguments are
therefore not interchangeable, and the pair of functions is what tells
them apart: `dbartsDrawLatents(fit = f, offset = o)` and
`dbartsDrawLatents(fit = f + o, offset = NULL)` draw from the same law
and agree, while `dbartsWorkingResponse` *differs* between the two, the
working response being stated in units the offset has been taken out of.
Pass the offset as `offset` to both.

**Why two functions.** For `"probit"`, `"ordinal"` and `"aft"` the
latent is a LOCATION and the working response is that location less the
offset. For `"logistic"` and `"nbinom"` it is a PRECISION - the
Polya-Gamma \\\omega_i\\ - and what a host regresses on is \\w_i (y_i -
1/2) / \omega_i - o_i\\ and \\(y_i - r) / (2 \omega_i) - o_i\\
respectively, under per-row precision weights \\\omega_i\\. For
`"student"` the latent is the scale-mixing \\\lambda_i\\, which weights
the row rather than entering the response, so the working response is
the gaussian one. A host that regresses on the latent itself is fitting
a different model, which is the composition defect these functions exist
to prevent.

**Random numbers.** The draws come from R's own stream, through a
generator created and destroyed within the call:
[`set.seed`](https://rdrr.io/r/base/Random.html) reproduces them, they
advance the stream, and they compose with any other R draw. They never
touch a sampler's chain generator, which is seeded separately and never
advances R's stream. That stream is *not* thread safe, so these
functions must not be called from a worker thread; call them on the main
thread and hand the result to the workers.

## Value

A numeric vector of length \\n\\. `dbartsDrawLatents` carries a
`"quantity"` attribute of `"location"` (`"probit"`, `"ordinal"`,
`"aft"`) or `"precision"` (`"logistic"`, `"nbinom"`, `"student"`),
matching what
[`getLatents`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
reports for that family.

## See also

[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md),
whose `getLatents`, `getFitsWithoutOffset` and `getDispersion` methods
report the same quantities from inside a running sampler.

## Examples

``` r
set.seed(99)
f <- rnorm(50L)
o <- rep(0.25, 50L)
y <- rbinom(50L, 1L, pnorm(f + o))

# one Albert-Chib step: the truncated-normal latent at the current fit, then
# the response a host regresses on
z <- dbartsDrawLatents("probit", fit = f, y = y, offset = o)
attr(z, "quantity")
#> [1] "location"
head(dbartsWorkingResponse("probit", latent = z, y = y, offset = o))
#> [1]  0.35539142 -0.62844846 -0.40419185  0.68244715 -0.03539326 -1.46411699

# logistic draws a PRECISION, so its working response is a quotient, not the
# latent shifted
w <- dbartsDrawLatents("logistic", fit = f, y = y, offset = o)
head(dbartsWorkingResponse("logistic", latent = w, y = y, offset = o))
#> [1]  3.238115 -4.067323 -7.838172  2.316351  4.837544 -2.638638
```
