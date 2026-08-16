# Embedding dbarts In A Larger Sampler

Reference for using BART as a conditional model inside a Gibbs or
Metropolis-Hastings scheme the caller writes: which channel an outer
block writes, which fit it reads back, whose prior is in force after a
swap, and how a composition is checked. The runnable set is
[`vignette("dbarts-as-a-component", package = "dbarts")`](https://vdorie.github.io/dbarts/articles/dbarts-as-a-component.md),
six named recipes with their output; this page is the rule each of them
applies.

## Details

A
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
is created once by
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) and then
run one sweep at a time with `run(0L, 1L)`, its inputs replaced between
sweeps. The fitting functions
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md),
[`bart2`](https://vdorie.github.io/dbarts/reference/bart.md) and
[`rbart_vi`](https://vdorie.github.io/dbarts/reference/rbart.md) own
their own loop and are not the surface for this; they return fit
objects, and the R modelling accessors
([`fitted`](https://rdrr.io/r/stats/fitted.values.html),
[`predict`](https://rdrr.io/r/stats/predict.html),
[`extract`](https://vdorie.github.io/dbarts/reference/bart.md)) apply to
those.

### The channels an outer block writes

Each method replaces one input in place, on every chain, taking effect
at the next sweep. All of them are documented in
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md).

|  |  |  |
|----|----|----|
| **method** | **what the outer block owns** | **typical use** |
| `$setResponse(y)` | the regression target | a partial residual, a working response, an imputed outcome |
| `$setOffset(o)` | an additive term BART does not fit | another block's linear predictor, a fixed exposure |
| `$setWeights(w)` | per-observation precision | a Polya-Gamma precision, a variance mixture, survey weights |
| `$setSigma(s)` | the residual standard deviation | a scale drawn elsewhere; requires `resid.prior = fixed()` |
| `$setPredictor(x, column)` | one predictor column | a latent covariate, an imputed predictor |
| `$setActiveRows(a)` | which rows enter the likelihood | a subset indicator drawn by another block |
| `$setCalibration(prior.scale = )` | the leaf prior in force | restating the prior after a response swap |

A predictor moved per observation rather than wholesale takes
`$setPredictor(x, column, forceUpdate = "partial")`, or
[`updatePredictorPerObservationJointly`](https://vdorie.github.io/dbarts/reference/updatePredictorPerObservationJointly.md)
when several samplers share the column. Both return a per-observation
install mask, `FALSE` where the new value would empty a leaf and was
rolled back. **That mask is part of a Metropolis accept decision, not a
diagnostic**: a declined observation holds its old value, so the
likelihood read after the call is its old likelihood, and treating the
resulting ratio as an acceptance leaves the host's copy of the latent
disagreeing with what the samplers hold.

Not every channel is open on every model. A refusal is part of the model
rather than a limitation - `setSigma` on a family that fixes the
residual scale by definition, a response swap that would re-anchor a
calibration stated at creation - and each names its reason.

### Which fit to read back

Five methods report a fitted quantity and no two report the same one;
the table in ‘Reading the fit’,
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md),
is authoritative. The rule an outer block needs: it conditions on
\\f(x_i)\\, so it reads `$getFitsWithoutOffset()` and adds back whatever
offset it installed. `run()$train` reports the same location with the
offset folded in, and `$getLatents()` reports the family's augmentation
variable, which carries no offset and for `"logistic"`, `"nbinom"` and a
Student-t residual distribution is a precision rather than a location at
all. Differencing the two is short by exactly the installed offset,
which reaches the outer block as a coefficient shrunk toward zero rather
than as an error.

### Whose prior is in force

The leaf prior is calibrated once, at creation, from the response the
sampler was built on. A per-sweep `setResponse` deliberately does *not*
re-anchor it (`updateScale = FALSE`, the default), so that the sweeps
are comparable - which also means a sampler built on a cold-start vector
keeps that vector's range as its prior scale for the whole run. State
the prior instead of inheriting it: `node.prior = normal(scale = )` at
creation (see
[`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)),
or `$setCalibration(prior.scale = )` afterward, with `$getCalibration()`
reporting what is actually in force.

A decomposition across \\K\\ samplers needs the same care in the other
direction. Each carries a prior sized for the whole of \\f\\, so their
sum has \\\sqrt{K}\\ times the prior standard deviation one sampler
would have; dividing each `prior.scale` by \\\sqrt{K}\\ restates them as
shares of one budget.

### The seeding contract

The value in force during a sweep is the one installed before it. Create
the sampler with the outer block's initial value already in place, and
write the channel before `run`, never after: a sweep run against a stale
offset, weight or response has drawn from the wrong conditional and
nothing later reports it. The mutators store no state unless
`updateState = TRUE` is passed explicitly, so a per-sweep write costs
nothing beyond the write itself.

### Augmentation written in R

The per-observation augmentation a non-gaussian sampler runs internally
is callable directly:
[`dbartsDrawLatents`](https://vdorie.github.io/dbarts/reference/dbartsAugmentation.md)
draws the latent at a supplied fit and
[`dbartsWorkingResponse`](https://vdorie.github.io/dbarts/reference/dbartsAugmentation.md)
turns it into the quantity a host regresses on. A host that needs the
latent for its own blocks runs the step itself over a gaussian sampler
with `resid.prior = fixed(1)`, rather than reimplementing the family.
Their `fit` argument is offset-free, as `$getFitsWithoutOffset()`
reports it.

### Checking a composition

Recovering a known truth on one simulated data set is weak evidence.
[`dbartsValidateComposition`](https://vdorie.github.io/dbarts/reference/dbartsValidateComposition.md)
runs simulation-based calibration over the composed one-sweep step and
reports, per scalar functional, whether the kernel targets the posterior
the stated prior and likelihood imply. It catches the defect class in
which a posterior mean is reported where a draw belongs.

## See also

[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)
for every method named here,
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) for
creating one,
[`dbartsAugmentation`](https://vdorie.github.io/dbarts/reference/dbartsAugmentation.md)
for the augmentation helpers,
[`dbartsValidateComposition`](https://vdorie.github.io/dbarts/reference/dbartsValidateComposition.md)
for the calibration check,
[`updatePredictorPerObservationJointly`](https://vdorie.github.io/dbarts/reference/updatePredictorPerObservationJointly.md)
for the shared-predictor sweep, and
[`samplePriorPredictive`](https://vdorie.github.io/dbarts/reference/samplePriorPredictive.md)
for inspecting a prior before fitting.
