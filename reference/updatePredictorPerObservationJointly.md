# Jointly Update a Shared Predictor per Observation Across Samplers

Applies a “partial” (per-observation) update of a single shared column
to several
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
objects at once, installing an observation's new value only when it
keeps every leaf node non-empty in every tree of every sampler.

## Usage

``` r
updatePredictorPerObservationJointly(samplers, x, column, updateState = NA)
```

## Arguments

- samplers:

  A
  [`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md),
  or a list of `dbartsSampler` objects, that share the same
  (index-aligned) observations. The samplers are treated as peers; none
  is privileged.

- x:

  A numeric vector of new values for the shared column, of length equal
  to the number of observations.

- column:

  A single integer or character string identifying the shared column; it
  cannot be missing. The column is matched *by name* across all
  samplers, so the shared variable may sit at different column positions
  in each. When given as an integer, it indexes the columns of the first
  sampler and the corresponding name is used to locate the column in the
  rest.

- updateState:

  A logical determining if the local cache of each sampler's state
  should be updated after the update completes. If `NA`, each sampler's
  own
  [`control`](https://vdorie.github.io/dbarts/reference/dbartsControl.md)
  default is used.

## Details

This is the multi-sampler analog of the `forceUpdate = "partial"` mode
of
[`setPredictor`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md).
A single sequential sweep installs each observation in every sampler at
once, and only if its new value keeps every leaf non-empty in every tree
of every chain of every sampler; otherwise that observation is rolled
back to its previous value in all of them. Like the single-sampler
partial mode, the update never changes tree structure – it only
re-routes observations.

This supports embedding several BART components that condition on a
common latent predictor (for example, a shared ability or trait in an
item-response model) within a larger Gibbs/Metropolis sampler.

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
