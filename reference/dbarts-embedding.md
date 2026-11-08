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
[`bartBT`](https://vdorie.github.io/dbarts/reference/bartBT.md) and
[`bart`](https://vdorie.github.io/dbarts/reference/bart.md) own their
own loop and are not the surface for this; they return fit objects, and
the R modeling accessors
([`fitted`](https://rdrr.io/r/stats/fitted.values.html),
[`predict`](https://rdrr.io/r/stats/predict.html),
[`extract`](https://vdorie.github.io/dbarts/reference/bartBT.md)) apply
to those.

### The channels an outer block writes

Each method replaces one input in place, on every chain, taking effect
at the next sweep. All of them are documented in
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md).

|  |  |  |
|----|----|----|
| **method** | **what the outer block owns** | **typical use** |
| `$setResponse(y)` | the regression target | a partial residual, a working response, an imputed outcome |
| `$setOffset(o)` | an additive term BART does not fit | another block's linear predictor, a fixed exposure |
| `$setWeights(w)` | per-observation precision, or under `weights.type = "frequency"` a count of identical observations | a Polya-Gamma precision, a variance mixture; under frequency, fractional counts such as expected memberships |
| `$setSigma(s)` | the residual standard deviation | a scale drawn elsewhere; use with `family = gaussian(sigma = fixed())` so the next `run()` does not redraw it |
| `$setPredictor(x, column)` | one predictor column | a latent covariate, an imputed predictor |
| `$setActiveRows(a)` | which rows enter the likelihood | a subset indicator drawn by another block |
| `$setLeafPrior(normal(sd = ))` | the leaf prior's spread | stating the prior a cold-start vector did not |

A predictor moved per observation rather than wholesale takes
`$setPredictor(x, column, forceUpdate = "partial")`, or
[`updatePredictorPerObservationJointly`](https://vdorie.github.io/dbarts/reference/updatePredictorPerObservationJointly.md)
when several samplers share the column. Both return a per-observation
install mask, `FALSE` where the new value would empty a leaf and was
rolled back; under a
[`monotone`](https://vdorie.github.io/dbarts/reference/monotone.md)
constraint the joint form also declines a value that would leave a
tree's leaf values out of order. **That mask is part of a Metropolis
accept decision, not a diagnostic**: a declined observation holds its
old value, so the likelihood read after the call is its old likelihood,
and treating the resulting ratio as an acceptance leaves the outer
sampler's copy of the latent disagreeing with what the samplers hold.

Not every channel is open on every model. A refusal is part of the model
rather than a limitation - `setSigma` on a family that fixes the
residual scale by definition, a response swap that would re-anchor a
calibration stated at creation - and each names its reason.

### When the augmentation variables are drawn

A sampler whose family augments its likelihood (the `probit` and
`ordinal` latents, the Polya-Gamma precisions of `logistic`, `binomial`
and `nbinom`, the Student-t scales, an `aft` sampler's imputed times)
draws those values at the start of each sweep, against the fit the
previous sweep left and the inputs in force, and then draws its trees
against them. A sweep therefore runs in this order: the augmentation
variables, then the trees, then sigma and the other parameters. No
setter draws them: `$setResponse`, `$setOffset`, `$setPredictor`,
`$setSigma`, `$setWeights`, `$setActiveRows`, `$setForestBasis`,
`$setData`, `$setModel` and `$setState` store what they are given, and
the next sweep draws the latents against the change.

Two kinds of outer step are both exact under this order. One conditions
on the latent: it runs the sampler, reads `$getLatents()`, and draws its
own quantity given those values. The other integrates the latent out,
scoring a proposal with the probit likelihood `pnorm(f)` of the
response, say, and then installing what it accepted; the next sweep
draws the latents under the accepted values, so the trees never fit a
latent drawn before the change. The order is also what makes a setter
call and a restore undo each other exactly: `$setWeights` with the old
weights and `$setState` with the stored state, in either order, leave
the sampler continuing bit for bit as one that installed the stored
state under the old weights would.

### Which fit to read back

Six methods report a fitted quantity and no two report the same one; the
table in ‘Reading the fit’,
[`dbartsSampler`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md),
is authoritative. The rule an outer block needs: it conditions on
\\f(x_i)\\, so it reads `$getFitsWithoutOffset()` and adds back whatever
offset it installed. `run()$train` reports the same location with the
offset folded in. The location latent of `"probit"`, `"ordinal"` and
`"aft"` that `$getLatents()` reports is on that same offset-inclusive
scale, \\z = f(x) + o + e\\, so the difference `z - train` is the
residual \\e\\, which carries none of the offset's signal. A block that
owns the offset regresses `z - $getFitsWithoutOffset()`; used the other
way, the offset reaches it as a coefficient shrunk toward zero rather
than as an error. The latent of `"logistic"`, `"nbinom"` and a Student-t
residual distribution is a precision, a location at no scale, with no
offset in it;
[`dbartsWorkingResponse`](https://vdorie.github.io/dbarts/reference/dbartsAugmentation.md)
gives the offset-free working response for any of them.

`$getLatents()` returns the values the last sweep's trees were drawn
from, so after a run it and `$getFitsWithoutOffset()` are a joint draw,
the latent being the value the fit was drawn from. Between a setter and
the next run it returns the values from before the setter (`$setData`
starts them as creation does, an ordinal or negative-binomial sampler
from the kept cutpoints or shape, and an `aft` `$setResponse` leaves a
censored row at its new bound), not values drawn against the change; the
next sweep draws them afresh.

### Whose prior is in force

The leaf prior is calibrated once, at creation, from the response the
sampler was built on. A per-sweep `setResponse` deliberately does *not*
re-anchor it (`updateScale = FALSE`, the default), so that the sweeps
are comparable - which also means a sampler built on a cold-start vector
keeps that vector's range as its prior scale for the whole run, and its
midrange as the prior's center (read back as
`$getLeafPrior()$prior.mean`, which no setter writes except through the
offset). State the prior instead of inheriting it:
`leaf.prior = normal(sd = )` at creation (see
[`dbartsPriors`](https://vdorie.github.io/dbarts/reference/dbartsPriors.md)),
or `$setLeafPrior(normal(sd = ))` afterward, with `$getLeafPrior()`
reporting what is actually in force. A caller whose rows enter and leave
the likelihood asks for the scale to be taken again from the rows then
in with `updateScale = TRUE` on `$setWeights` or `$setActiveRows`, in
burn-in only, and a `linear` or `gp` leaf's covariate center and scale
with `updateStandardization = TRUE` on the same call; `FALSE` after
burn-in holds both, so every kept draw comes from one prior (see
[`updateStandardization`](https://vdorie.github.io/dbarts/reference/dbartsSampler-class.md)).
A named sd is absolute: a channel that does re-anchor the transform
leaves it where it was.

A decomposition across \\K\\ samplers needs the same care in the other
direction. Each carries a prior sized for the whole of \\f\\, so their
sum has \\\sqrt{K}\\ times the prior standard deviation one sampler
would have; dividing each `prior.sd` by \\\sqrt{K}\\ restates them as
shares of one budget.

### The seeding contract

The value in force during a sweep is the one installed before it. Create
the sampler with the outer block's initial value already in place, and
write the channel before `run`, never after: a sweep run against a stale
offset, weight or response has drawn from the wrong conditional and
nothing later reports it. The mutators store the state after each write
when `control@updateState` is `TRUE`, the default, so a sampler driven a
sweep at a time is built with `dbartsControl(updateState = FALSE)`, as
the vignette's recipes are; a per-sweep write then costs nothing beyond
the write itself.

### Augmentation written in R

The per-observation augmentation a non-gaussian sampler runs internally
is callable directly:
[`dbartsDrawLatents`](https://vdorie.github.io/dbarts/reference/dbartsAugmentation.md)
draws the latent at a supplied fit and
[`dbartsWorkingResponse`](https://vdorie.github.io/dbarts/reference/dbartsAugmentation.md)
turns it into the quantity the outer sampler regresses on. An outer
sampler that needs the latent for its own blocks runs the step itself
over a gaussian sampler with `family = gaussian(sigma = fixed(1))`,
rather than reimplementing the family. Their `fit` argument is
offset-free, as `$getFitsWithoutOffset()` reports it.

### Driving it from compiled code

A package that runs its outer block in C or C++ declares
`LinkingTo: dbarts` and includes `dbarts/dbarts.h`, the flat C API. It
builds the sampler in R from its own R code and then reads the handle
out of that object's external pointer with `R_ExternalPtrAddr`
(`$getPointer()` is what hands the pointer down).
[`dbarts`](https://vdorie.github.io/dbarts/reference/dbarts.md) takes a
formula or a matrix, not a model object, so a consumer starting from its
own already-built `(control, model, data)` triple
([`dbartsSpec`](https://vdorie.github.io/dbarts/reference/dbartsSpec.md)
resolves one without constructing a sampler) builds the sampler instead
with `methods::new("dbartsSampler", control, model, data)`, then calls
`$getPointer()` on the result. There is no creation entry point and no R
type in the header; the handle is the whole boundary.

The handle is valid until the R object is garbage collected or replaces
its pointer, which it does whenever it re-creates its engine from a
stored state - after a [`load`](https://rdrr.io/r/base/load.html), or
after `$setState`. Re-read the handle after any such restore, and keep
the R object reachable for as long as the C side holds one.

What the header reaches is the per-sweep loop: `dbarts_sampler_run` into
caller-owned buffers, the conditioning channels `setResponse`,
`setOffset` and `setSigma`, `getLatents`, `predict`, tree storage and
printing, and the queries the outer sampler sizes its results struct
from, among them `dbarts_sampler_numFittedValuesPerObservation` and
`dbarts_sampler_numVariableCountForests`, the widths of a multinomial or
multi-forest sampler's fits and split counts; a flat run polls for R's
interrupt as `$run` does and, once the workers have joined, raises it as
R's interrupt condition (not an error:
[`try()`](https://rdrr.io/r/base/try.html) lets it through), leaving the
sampler usable, with any saved results of that call to be discarded and
a pass over the trees possibly partly applied; it also warns of a slow
monotone leaf-order count once per sampler. Everything else in the table
above - the predictor and weight channels, the active-row mask, the
state round trip, tree extraction, the multi-forest surface - is an R
method on the same object, and the two views drive one engine. The
setters COPY what they are handed into buffers the sampler owns at
creation, so the caller's array is free the moment the call returns and
conditioning on new values means calling the setter again rather than
writing through the array.

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
