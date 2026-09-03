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
every sufficient statistic and every family-level parameter update, and keeps
its leaf occupancy and its fitted value ([[R/dbarts.R#setActiveRows]],
[[man/dbartsSampler-class.Rd#setActiveRows]]). Whether its own latent is drawn
as well is a per-family matter, settled under "The mechanism".

Zero case weights already do this for a gaussian response - a masked gaussian
sampler is bitwise `setWeights(w * a)` with no mask - and for Student-t through
the same composite. What the channel adds is the families zero weights cannot
reach. Probit, ordinal, negative binomial and AFT carry no case weights at all:
each refuses them at creation, and refuses a post-creation change by name
([[src/R_interface_bartcore.cpp#refuseBinaryWeightChange]]), each for its own
reason: a weighted truncated-normal latent likelihood is not a coherent model
for probit and ordinal, count exposure belongs in the offset for negative
binomial, and AFT is a gaussian on log-time that simply carries none. Logistic
does take a weight change, but its weights are Polya-Gamma copy counts, which
must be positive integers on every surface, so no count expresses "not in the
data set" ([[src/R_interface_bartcore.cpp#enforceBinaryWeightPolicy]]).
Multinomial's response is a count matrix with no weight channel at all.

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

The mask and the case weights are absolute and independent: installing one
never disturbs the other, and a family that holds both serves their composition
whichever arrived first. What that composition is varies by family - `w * a`
for gaussian, `w * lambda * a` for Student-t, the mask alone for probit and
ordinal, `a * omega` for logistic and negative binomial - and "Per family"
below gives each. `setResponse` and `setOffset` leave an installed mask
standing, because it names rows and not responses; `setData` clears it, because
the mask is length-n and n may change.

## The mechanism

One validating and normalizing scan owns the channel, in the engine rather than
in any host, so every surface inherits it
([[src/bartcore/chain.hpp#Chain::setActiveRows]]). It refuses when the family
implements no mask, scans for an element that is neither 0 nor 1, and
normalizes an all-ones mask to a null pointer before any family sees it. The
values are copied; nothing downstream retains the caller's array. That first
refusal is family-generic and no shipped family reaches it, multinomial
included; it can only be raised by a family added later that does not override
the base
([[src/R_interface_bartcore.cpp#"active-row masking is not implemented for this response family"]]).

Below that scan each family composes `a` into *its own* precision vector,
rather than the chain holding a separate mask every kernel would consult
([[src/bartcore/model.hpp#ResponseModel::setActiveRows]]). A zeroed precision
is what drops the row from every leaf sufficient statistic, branch score and
leaf draw, with no edit to the tree machinery.

Whether the mask also suppresses a family's latent draw is a per-family
question, settled by the primitive. Where the augmentation is drawn by
rejection - the truncated normal behind probit, ordinal and AFT, the
Polya-Gamma variate behind logistic, negative binomial and multinomial - the
number of uniforms consumed depends on the argument, so an inactive row's draw
is skipped rather than taken and discarded: a discard would desynchronize the
chain's generator from a sampler built on the retained rows alone. That row's
latent keeps its last drawn value, which stays finite, and is stale until the
row is active again
([[src/bartcore/model.hpp#ProbitResponse::refreshLatents]]). Student-t is the
deliberate exception. Its lambda comes from a gamma whose consumption depends
on the shape `(nu + 1) / 2` and not on the mask, so lambda is drawn at every
row and the mask annihilates it through the composite instead, at no cost to
the generator ([[src/bartcore/model.hpp#TResponse::refreshLatents]]).

## Per family

Gaussian multiplies the case weight, so the sigma posterior's degrees of
freedom recount off the composite rather than staying at the unmasked total
([[src/bartcore/model.hpp#GaussianResponse::setActiveRows]]).

Student-t multiplies the mask into the same scale-mixture composite,
`c_i = w_i lambda_i a_i`, and hands that to its contained gaussian, so the node
statistics, the sigma degrees of freedom and the nu grid statistics all inherit
it - an inactive row leaves them exactly as a zero user weight would
([[src/bartcore/model.hpp#TResponse::setActiveRows]]).

Probit carries no user weights, so the mask *is* the precision vector
([[src/bartcore/model.hpp#ProbitResponse::workingWeights]]). Ordinal serves it
the same way, and recomputes the cutpoint proposal scales on install because
they are derived from the category counts, which the mask changes
([[src/bartcore/model.hpp#OrdinalResponse::setActiveRows]]). The cutpoint
acceptance target reads the mask live and needs no install-time work.

Logistic and negative binomial serve a *separate* composite `a * omega` rather
than writing the zero into omega itself: the working response divides by omega,
and 0 times infinity in the node kernels is a NaN
([[src/bartcore/model.hpp#LogisticResponse::workingWeights]]). Negative
binomial also restricts the collapsed statistic the dispersion grid draw reads,
and rebuilds the count histogram the grid's own kernel is built from over the
active rows at every mask change - the channel's only per-install cost above a
pass over the rows ([[src/bartcore/model.hpp#NBResponse::setActiveRows]],
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

No decorator needs a mask of its own, because each reads a precision the
response has already composed. An additive coupling takes the inert default
([[src/bartcore/combiner.hpp#ForestCombiner::setActiveRows]]) and picks the mask
up from the precision the sweep hands its per-forest response, which scales it
by the squared forest multiplier
([[src/bartcore/combiner.hpp#AmplitudeForestCombiner::formForestResponse]]).
The heteroscedastic mean weights divide that same composed precision by the
variance surface ([[src/bartcore/chain.hpp#formMeanWeights]]). A grouped
random-intercept response delegates to its base family
([[src/bartcore/model.hpp#GroupedResponse::setActiveRows]]).

## Reporting, and what the mask is not

The pointwise log-likelihood reports NaN at an inactive row - the channel's own
"not in the model" flag - rather than the finite value the row's fit would
still give ([[src/bartcore/model.hpp#GaussianResponse::computeLogLikelihood]]),
so anything summing that channel must decide what to do with the NaNs.

The mask is not saved state. The serialized state carries derived quantities
and no raw conditioning vector, and the reference-class sampler mirrors no mask
onto its data object, so a sampler RE-CREATED from a stored state - what
`$getPointer` does when the external pointer has gone stale across a save and
load - has none and must be given one again
([[docs/design/bart-as-a-component.md#What engine state does not carry]]).
Installing a state into a live sampler is a different thing:
[[src/bartcore/chain.hpp#Chain::setState]] touches no mask, so one already
installed stays in force across `$setState`. Nor is there a getter, so what is
installed is the caller's own record.

## The surfaces

Three, all reaching the same scan. The reference-class method
[[R/dbarts.R#setActiveRows]] validates length and values in R first, for the
message. The bridge entry
[[src/R_interface_bartcore.cpp#bartcore_setActiveRows]] and the flat C entry
[[src/C_interface.cpp#dbarts_sampler_setActiveRows]] each probe
[[src/bartcore/facade.hpp#SamplerShape::supportsActiveRows]] and never the
family, since a Student-t sampler reports as gaussian and does accept a mask.
That probe is read off the family, chain by chain, from the same predicate the
setter refuses on
([[src/bartcore/chain.hpp#Chain::supportsActiveRows]]), so the advertised
capability and the refusal cannot disagree, and
[[src/bartcore/sampler.hpp#Sampler::setActiveRows]] fans out to every chain
under a refusal they all share, so an install cannot land half applied.

## The tests that pin it

Per family, the masked n-row kernel against the same kernel run over the
compacted active rows, bitwise in value and in variates consumed - which is
what tells a skipped draw from a discarded one
([[tests/cpp/test_model.cpp#testActiveRowsProbitKernel]],
[[tests/cpp/test_model.cpp#testActiveRowsOrdinalKernels]],
[[tests/cpp/test_model.cpp#testActiveRowsLogisticKernel]],
[[tests/cpp/test_model.cpp#testActiveRowsNBKernels]],
[[tests/cpp/test_model.cpp#testActiveRowsAFTCensored]],
[[tests/cpp/test_sampler.cpp#testActiveRowsMultinomialKernel]]). Gaussian and
Student-t are pinned instead against a bare arm carrying `w * a` as fixed
weights, the degrees-of-freedom recount included
([[tests/cpp/test_model.cpp#testActiveRowsGaussianDf]],
[[tests/cpp/test_model.cpp#testActiveRowsStudentComposite]]).

At sampler level: the normalizer, the value refusal, the all-zeros run and a
gaussian arm against `setWeights(w * a)`
([[tests/cpp/test_sampler.cpp#testActiveRows]]); the same on a forest already
grown into the rows the mask removes
([[tests/cpp/test_sampler.cpp#testActiveRowsOnGrownForest]]); substituted
responses at the inactive rows leaving every active row's recorded draw bitwise
([[inst/tinytest/test-active-rows-pins.R#"logistic, nbinom and aft"]],
[[inst/tinytest/test-active-rows-pins.R#"multinomial, GLOBAL only"]]); the
gaussian, Student-t, BCF, heteroscedastic and grouped arms bitwise against
`setWeights(w * a)` ([[inst/tinytest/test-active-rows-pins.R#"heteroSampler"]]);
and the flat C pins for the all-ones no-op, the NULL clear and the fractional
refusal ([[inst/tinytest/test-capi.R#"capi_set_active_rows"]]).
