# The active-row mask

Status: reference, current. Documents the shipped channel as-is; not a design
proposal and carries no landing date - update in place whenever the channel
changes.

A sampler holds a per-observation 0/1 vector saying which rows are in its data
set this sweep.

## What it is for

An outer sampler that redraws a row subset between sweeps - principal
stratification, mediation, an instrumental-variable compliance class - has to
change which rows a BART block conditions on without rebuilding it. Rebuilding
throws away the tree structure, the latents and the generator position;
compacting to the retained rows throws away every dropped row's fitted value.

`$setActiveRows(a)` changes membership in place. A row with `a[i] = 0` leaves
every sufficient statistic, every family-level parameter update and its own
latent draw, and keeps its leaf occupancy and its fitted value
([[R/dbarts.R#setActiveRows]], [[man/dbartsSampler-class.Rd#setActiveRows]]).

Zero case weights already do this for a gaussian response, and a masked
gaussian sampler is bitwise `setWeights(w * a)` with no mask. What the channel
adds is the latent families, where a weighted likelihood is not a coherent
model and a post-creation weight change is refused outright: probit, ordinal,
logistic, negative binomial, AFT and multinomial.

## The contract

`active` is a length-n vector of doubles, each element exactly 0 or 1. A null
pointer, or `NULL` from R, clears whatever is installed.

An all-ones mask is accepted and installs *nothing*: the family restores its
pre-mask precision pointer by identity, so a mask returning to all ones is the
same object as no mask at all. `setWeights` is the opposite: an all-ones vector
there does install. One channel is membership, the other precision.

An all-zeros mask is accepted and runs. Every forest sits at its prior and
every row still receives a fit, so a stratum that empties needs no special case
in the caller.

A fractional element is refused rather than rounded or honoured: a value
strictly between 0 and 1 is a weighted likelihood, which the latent families
have no coherent form for, and belongs to `setWeights`
([[src/R_interface_bartcore.cpp#refuseNonBinaryMask]],
[[inst/include/dbarts/dbarts.h#dbarts_sampler_setActiveRows]]). NaN fails both
equality tests, so it is refused too.

The mask and the case weights are absolute and independent: a family holds both
and serves `w * a` in either call order. `setResponse` and `setOffset` leave an
installed mask standing, because it names rows and not responses; `setData`
clears it, because the mask is length-n and n may change.

## The mechanism

One validating and normalizing scan owns the channel, in the engine rather than
in any host, so every surface inherits it
([[src/bartcore/chain.hpp#Chain::setActiveRows]]). It refuses when the family
implements no mask, scans for an element that is neither 0 nor 1, and
normalizes an all-ones mask to a null pointer before any family sees it. The
values are copied; nothing downstream retains the caller's array.

Below that scan each family composes `a` into *its own* precision vector,
rather than the chain holding a separate mask every kernel would consult
([[src/bartcore/model.hpp#ResponseModel::setActiveRows]]). A zeroed precision
is what drops the row from every leaf sufficient statistic, branch score and
leaf draw, with no edit to the tree machinery.

The latent draw at an inactive row is skipped, not drawn and discarded. These
augmentations - the truncated normal, the Polya-Gamma variate - are rejection
samplers whose consumption depends on their argument, so drawing and discarding
would desynchronize the stream from a sampler built on the retained rows. The
row's latent keeps its last drawn value, which stays finite, and is stale until
the row is active again
([[src/bartcore/model.hpp#ProbitResponse::refreshLatents]]).

## Per family

Gaussian multiplies the case weight, so the sigma posterior's degrees of
freedom recount off the composite rather than staying at the unmasked total
([[src/bartcore/model.hpp#GaussianResponse::setActiveRows]]).

Probit carries no user weights, so the mask *is* the precision vector
([[src/bartcore/model.hpp#ProbitResponse::workingWeights]]). Ordinal serves it
the same way and additionally recomputes the cutpoint proposal scale on
install, since both cutpoint sums and the acceptance target restrict to the
active rows ([[src/bartcore/model.hpp#OrdinalResponse::setActiveRows]]).

Logistic and negative binomial serve a *separate* composite `a * omega` rather
than writing the zero into omega itself: the working response divides by omega,
and 0 times infinity in the node kernels is a NaN
([[src/bartcore/model.hpp#LogisticResponse::workingWeights]]). Negative
binomial also restricts the collapsed statistic the dispersion grid draw reads
and rebuilds the count-histogram kernel behind it over the active rows at every
mask change, which is the channel's one per-install cost
([[src/bartcore/model.hpp#NBResponse::setActiveRows]],
[[src/bartcore/model.hpp#NBDispersionPrior::computeKernel]]).

AFT composes into its contained gaussian, inheriting that recount, and skips an
inactive censored row's log-time redraw; the response transform stays the
full-data one ([[src/bartcore/model.hpp#AFTResponse::setActiveRows]]).

Multinomial's mask is global, and lands on the softmax coupling rather than the
response, which holds no precisions of its own
([[src/bartcore/model.hpp#MultinomialResponse::supportsActiveRows]],
[[src/bartcore/combiner.hpp#MultinomialForestCombiner::setActiveRows]]). An
inactive row's K interleaved Polya-Gamma draws are skipped
([[src/bartcore/combiner.hpp#drawForestGlue]]) and its composed precision is
zeroed in every category
([[src/bartcore/combiner.hpp#MultinomialForestCombiner::formForestResponse]]),
while the row keeps its leaf occupancy and its reported probabilities. A
per-forest mask is refused permanently, on model grounds: the softmax margin is
a log-sum-exp over the other K-1 forests, so "row i is out of category f only"
restricts no likelihood at all.

The decorators carry no mask code of their own. An additive coupling takes the
inert default ([[src/bartcore/combiner.hpp#ForestCombiner::setActiveRows]]) and
inherits the composed precision
([[src/bartcore/chain.hpp#composeForestWeights]]); the heteroscedastic mean
weights divide that composite by the variance surface
([[src/bartcore/chain.hpp#formMeanWeights]]); a grouped random-intercept
response delegates to its base family
([[src/bartcore/model.hpp#GroupedResponse::setActiveRows]]).

## What it does not touch

The pointwise log-likelihood reports NaN at an inactive row - the channel's own
"not in the model" flag - rather than the finite value the row's fit would
still give ([[src/bartcore/model.hpp#GaussianResponse::computeLogLikelihood]]),
so anything summing that channel must decide what to do with the NaNs.

The mask is not saved state. The serialized state carries derived quantities
and no raw conditioning vector, and the R5 sampler mirrors no mask onto its
data object, so a re-created or state-restored sampler has none and must be
given one again
([[docs/design/bart-as-a-component.md#What engine state does not carry]]).

## The surfaces

Three, all reaching the same scan. The reference-class method
[[R/dbarts.R#setActiveRows]] validates length and values in R first, for the
message. The bridge entry
[[src/R_interface_bartcore.cpp#bartcore_setActiveRows]] and the flat C entry
[[src/C_interface.cpp#dbarts_sampler_setActiveRows]] each probe
[[src/bartcore/facade.hpp#SamplerShape::supportsActiveRows]] and never the
family, since a Student-t sampler reports as gaussian and does accept a mask.
That probe is derived from the predicate the setter refuses on
([[src/bartcore/facade.hpp#SamplerBase::setActiveRows]]), so the advertised
capability and the refusal cannot disagree, and
[[src/bartcore/sampler.hpp#Sampler::setActiveRows]] fans out to every chain
under a refusal they all share, so an install cannot land half applied.

## The tests that pin it

Per-family kernel comparisons against the compacted case, bitwise in value and
in generator stream ([[tests/cpp/test_model.cpp#testActiveRowsNBKernels]],
[[tests/cpp/test_sampler.cpp#testActiveRowsMultinomialKernel]],
[[tests/cpp/test_sampler.cpp#testActiveRows]]). At sampler level, substituting
arbitrary responses at the inactive rows and requiring every active row's draw
to stay bitwise, which fails outright if a skipped draw is taken and discarded
instead ([[inst/tinytest/test-active-rows-pins.R#"logistic, nbinom and aft"]],
[[inst/tinytest/test-active-rows-pins.R#"multinomial, GLOBAL only"]]). And the
flat-C pins for the all-ones no-op, the NULL clear and the fractional refusal
([[inst/tinytest/test-capi.R#"capi_set_active_rows"]]).
