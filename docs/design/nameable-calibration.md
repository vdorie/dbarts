# Naming the leaf-prior calibration

Status: PARTIAL, 2026-08-13. Creation half LANDED c2a7e89b, mid-chain half
LANDED (this commit); the flat-C half is an item inside the dbarts.h reshape's
S1. Plan: docs/plans/nameable-calibration.md. Companion:
docs/design/prior-defaults.md, which owns the defaults this names an
alternative to, and docs/design/r-c-division.md, whose "the wall that remains
is calibration" is what this closes.

An R program can state the leaf prior it composes against, in response units,
instead of inheriting one from the range of whatever vector the sampler
happened to be constructed on. A structurally correct pure-R probit
composition sweeping its construction range 16x moves the implied leaf sd 16x
and the posterior sd of `f` by about 4.6x, with no error and no warning; the
composed model's blocking ports to R, but its prior does not unless it can be
named.

## 1. The named quantity

For a forest with `m` trees, internal leaf scale `s`, precision divisor `k`,
and response transform `(fitScale, fitShift)`:

    prior.scale = fitScale * s * sqrt(m)    # response units, forest-total sd at k = 1
    prior.sd    = prior.scale / k           # response units, at the current k
    prior.mean  = fitShift                  # response units

`prior.scale` is the identified object. Only the ratio `s / k` enters any draw
law, so under a fixed `k` the pair `(prior.scale, k)` has one degree of
freedom; `k` matters when and only when a hyperprior draws it, and then
`prior.scale` is the constant of the model while `prior.sd` moves every sweep.
`prior.mean` is a read-only diagnostic whose lever is the offset channel.

Three priors are deliberately NOT part of this: the residual-variance prior
(its anchor is the named `dbarts(..., sigma =)`, already original-scale end to
end), the tree prior (scale free), and the variance forest's own leaf prior
(an inverse-gamma-shaped `ConstantVarianceLeaf` that sits outside `forests_`
and is not addressable - a door, not coverage).

## 2. The surface

Creation, `docs/design/prior-defaults.md` for the defaults it overrides: a
`prior.scale` slot on `dbartsModel`, spelled `node.prior = normal(scale =)` -
also `linear(scale =)` and `gp(scale =)`, and `sd =` for the same quantity at
the resolved `k` - and a `prior.scale` formal of `bart`, `bart2` and
`rbart_vi`. Unset it is `NA` and nothing changes.

Mid-chain, over every chain:

    cal <- sampler$getCalibration(forest = 1L)
    # numeric matrix, one ROW per chain, columns prior.scale prior.sd
    # prior.mean k k.has.hyperprior response.scale response.shift
    # amplitude.prior.variance amplitude.prior.scale node.scale.factor
    # node.scale.divisor basis.row.norm,
    # plus attr(, "leaf.model") in {constant, monotone, linear, gp}

    sampler$setCalibration(prior.scale = 3.0, forest = 1L, updateState = NA)
    sampler$setCalibration(prior.sd = 1.5, forest = 1L)   # sugar; fixed k only

The engine side is two entries on `Chain` - `forestCalibration` and
`setForestPriorScale` - sharing one private `priorScaleFactor` helper, so the
reader and the writer are the same conversion in both directions rather than
two spellings of it. Both are total over the four leaf models (each carries
the one `scale` field) and carry no family switch: the transform is one
virtual call, so the per-family units below are a documentation table, not a
branch.

| family | units of prior.scale / prior.mean | response.scale | response.shift |
|---|---|---|---|
| gaussian | response | `range_` | `range_*0.5 + min_` |
| aft | log survival time, anchored on the OBSERVED log times at creation and never re-anchored as censored latents are imputed | `range_(log T)` | `range_*0.5 + min_` |
| Student-t | response (delegates) | `range_` | as gaussian |
| grouped | response (delegates); the tau prior is a separate object | `range_` | as base |
| probit, weighted binary | probit latent | 1 | 0 |
| ordinal | probit latent, relative to the pinned first cutpoint | 1 | 0 |
| logistic | logistic latent (log-odds) | 1 | 0 |
| nbinom | log-odds `psi`; a mean-scale reading needs the current dispersion | 1 | 0 |
| multinomial | softmax log-odds, per category forest | 1 | 0 |
| BCF | response; `k` fixed at 1 by the map | `range_` | `range_*0.5 + min_` |
| BCF, probit | probit latent; `k` fixed at 1 by the map | 1 | 0 |
| BCF, logistic | logistic latent (log-odds); `k` fixed at 1 by the map | 1 | 0 |

On every BCF row, `prior.scale` is the prior sd of forest f's own `f_f`, NOT of
its contribution `m_f f_f` to the combined index - `$getForestAmplitudes()`
reports the multiplier, the other half of the pair. It is stated PER UNIT OF
BASIS ROW NORM; under the two latent rows the index itself is in latent sd
units with sigma PINNED, where the gaussian row's is `sd(y)` and a drawn sigma
partly absorbs a mis-scaled basis. The induced index prior, the map's `sqrt(K)`
dispersion and the `sqrt(2/K)` node scale factor default that cancels it are
docs/design/multiplier-combiner.md's map section. Two of the five columns are
therefore where a defaulting caller READS the shipped prior: after
`binary-kforest-prior-default` S2, `node.scale.factor` prints `sqrt(2/K)` and
`amplitude.prior.scale` prints the family's own half-Cauchy median, 2 under
gaussian and 1 under the latent families.

The five map columns are the same rows read one level down
(`binary-kforest-prior-default` S1). On a mapped forest they carry the
DECOMPOSITION of `prior.scale`, `factor * s / (divisor * rowNorm)`, so the
anchor s - the one map quantity with no column, and data-dependent under
gaussian - is recovered as `prior.scale * divisor * rowNorm / factor`; every
other row above reports NaN in all five, which is a positive "not map-derived"
signal rather than a neutral 1.0. Two of them are the amplitude PRIOR, in
whichever of `ForestAmplitudePrior`'s two exclusive spellings the forest
carries, and never the scale mixture's live variance auxiliary, which is state
and rides `$state`'s `aVariance`. The engine reads its own retained spec, not
the combiner, so no virtual joined `ForestCombiner` and no draw moved.

| channel | the five map columns |
|---|---|
| `setForestBasis` | `basis.row.norm` re-derived from the new basis, the other four unchanged, and the two `node.scale` columns RESTORED if an install had cleared them |
| `setState` / `installTrees` with a foreign leaf scale | `node.scale.factor` and `node.scale.divisor` go NaN (the stored pair no longer decomposes what is in force); `amplitude.prior.variance` FOLLOWS the state; `basis.row.norm` unchanged |
| `setState` restoring a sampler's OWN state | unchanged, since the installed scale is bitwise the one in force |

## 3. What prior.sd means, per leaf model

> `prior.scale` and `prior.sd` describe the LEAF-PARAMETER scale of the forest
> total. They equal the prior sd of `f(x)` at every x for the constant leaf
> only. For the other three the prior of `f(x)` is x-dependent and `prior.sd`
> bounds it in a leaf-model-specific direction.

- **constant**: EXACT, and `prior.mean` is the prior mean of `f(x)`.
- **linear**: a LOWER bound, attained at the standardized covariate origin;
  `sd(f(x)) = prior.sd * sqrt(1 + ||z(x)||^2)`. `prior.mean` exact.
- **gp**: an UPPER bound over x, attained at rows reproducing a leaf member
  and on over-cap leaves; elsewhere `prior.sd^2 c(x)' C^-1 c(x)`, decaying to
  0 as x leaves the leaf's cloud, where every draw equals `prior.mean`.
- **monotone**: a LOWER bound in the interior (realized sd a few per cent to
  ~20% above it), and `prior.mean` is NOT the prior mean of `f(x)` under an
  active constraint - that marginal is skew, with an x-dependent mean tracking
  the constraint direction across several `prior.sd`.

The setter is TOTAL over all four - it writes the one `scale` field each
carries, and that write is exactly well defined in every case - while the read
is bounded for three of them. The asymmetry is stated rather than hidden, and
`inst/tinytest/test-calibration-prior-draws.R` measures each bound in the
direction it claims.

## 4. Authority

`model@prior.scale` records the NAMED INTENT, in response units. The ENGINE
holds what is IN FORCE, and `$getCalibration` is the only authoritative
reader. The intent is applied at creation and at every `setModel`, and is
never rewritten by the engine. A channel that re-anchors the transform
(`setResponse` / `setOffset` at `updateScale = TRUE`, `setData`) moves what is
in force without touching the intent, and the getter shows the move -
re-anchoring stops being silent, which is the whole point of shipping a reader
beside the writer. Re-issuing the intent is `$setCalibration` or
`$setModel(sampler$model)`.

Per channel:

| channel | prior.scale in force |
|---|---|
| `setResponse` / `setOffset` at `updateScale = FALSE`, `setWeights`, `setSigma` | unchanged |
| `setResponse` / `setOffset` at `updateScale = TRUE`, `setData` | re-anchored, so it moves |
| `setModel` | re-derived from the model's `prior.scale` against the CURRENT transform when finite; otherwise replaced through `node.scale`, internal units. Also re-pins sigma for gaussian/aft with no variance forest |
| `$setCalibration` | written, every chain; nothing else moves |
| `storeState` / `setState` | adopted from the state, which carries the leaf scale |
| warm start (`installTrees`) | ADOPTED from the donor, whose trees were drawn under its scale. Recipe: re-issue `$setCalibration` afterwards to keep your own |

## 5. Refusals, and why each is one

1. `prior.mean =` is refused: the leaf values it would shift are already
   drawn. The lever is `setOffset(rep_len(-prior.mean, n))`.
2. `prior.sd =` under a sampled `k` is refused, naming both remedies (write
   `prior.scale`, or pin `k` through `setModel`). The binary families default
   to a `k` hyperprior, and a named sd under one drifts silently. It is also
   refused when the chains' `k` have diverged, since one number would name a
   different scale on each; `prior.scale` is k-free and serves both cases.
3. BCF and multinomial forests are refused, at three creation sites and again
   mid-chain: their per-forest leaf scales come from their own calibration
   maps, which would drop a named value in silence rather than honor it.
4. Host-shell samplers are refused through `refuseHostMutation`.
5. A non-finite or non-positive value is an ERROR, not a refusal. `NaN` is
   among them even though `is.na(NaN)` is TRUE in R: it names no intent and
   cannot serve as a divisor.
6. DART samplers are NOT refused, which is the concrete gain over the existing
   route - `$setModel` refuses a DART sampler outright, while this write
   touches no split machinery.
7. The variance forest is not addressable; it is a separate leaf model outside
   `forests_`. A door, recorded.

## 6. Exactness

The reader is a pure function of engine state, so the writer must not perturb
it. `setForestPriorScale` SKIPS the write when the requested value reproduces
what is in force, on either spelling of it - the internal scale it derives,
and the `prior.scale` the reader reports - so a read-then-write is bitwise
inert and cannot move the equivalence baseline. That is the load-bearing
direction, and it is measured: `inst/tinytest/test-calibration-midchain.R`
runs identical sweeps across a get-then-set, on a named scale and on an
inherited one.

The other direction is pinned at ulp level rather than bitwise, and that is a
robustness choice rather than a concession to a failure the fixture shows. The
value in force is the request divided by the transform, and the report
multiplies it back, so exactness is a property of the particular `(P, f)` pair:
`(P / f) * f != P` for about 10% of positive pairs (MEASURED 194596 of 2e6).
The shipped fixture is mostly on the lucky side of that - at its own transform
the round trip is BITWISE EXACT for all four requests at `n.trees` 20, 50 and
200, and the one rounding cell among them is `m = 25` with `P = 0.25`, off by
half an ulp - which is exactly why the assertion is not written as bitwise: it
would be pinning an accident of this response vector and this tree count, and
any fixture edit, any re-anchoring channel, or any consumer's own transform
could break it without a defect. Recovering the bit in general would mean
either caching the named value in the engine - which the state format and
every re-anchoring channel would then have to carry, and which would make the
reader something other than a reader - or nudging the scale in force off the
exact quotient. Neither is worth an ulp, so set-then-get fidelity is pinned at
ulp level, fifteen orders of magnitude below the `sqrt(m)`-shaped error it
exists to catch, and the O4b loop carries the rounding cell so the tolerance
is exercised rather than merely permitted.
