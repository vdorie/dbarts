# Forest combiner: design

Status: LANDED, 2026-07-14. Promotes the BCF glue that lived as a Chain-side
special case (bcf_, drawGlue, combinedFits, formForestResponse,
forestMultiplier, the if(bcf_) sweep branches - the shape bcf.md's Forest
split shipped in) into a polymorphic ForestCombiner<L> hierarchy beside
Forest<L> (src/bartcore/combiner.hpp since the multinomial extraction below;
src/bartcore/chain.hpp at this refactor's own landing). BCFForestCombiner<L> is
its first instance, the math it carried unchanged from bcf.md and
bcf-ridge-interweaving's landing at the time - since GENERALIZED to the K-forest
basis/amplitude family, whose math is docs/design/multiplier-combiner.md's;
MultinomialForestCombiner<L>
(docs/design/multinomial.md) is now the second. docs/plans/forest-combiner.md carries the
step plan, its binding contracts, and its resolved Open questions; this note
records the shape as landed and what it does and does not anticipate.

## Goal

A Chain that holds more than one forest (BCF today) needs somewhere to put
per-forest residual formation, the coupling draw, the combined per-observation
location, and per-forest reporting addressing. Before this refactor that
somewhere was Chain itself, hardcoded to forest 0 = mu / forest 1 = tau at
every touchpoint - workable for exactly one multi-forest model, but a wall for
the next one (multinomial, heteroscedastic, hurdle) and for composing a
multi-forest model with GroupedResponse. The combiner gives multi-forest
coupling an owner with a virtual surface a second model can implement, while
costing the single-forest chain nothing: combiner_ is null off BCF, and every
touchpoint collapses to the direct forest-0 path with no virtual call.

## The ForestCombiner<L> hierarchy

ForestCombiner<L> (chain.hpp, beside Forest<L>) is a virtual base templated on
the leaf model, for the same reason Forest<L> is: the post-combine move reaches
Forest<L>'s tree/fit buffers and saved-tree FlatNodes directly. A Chain holds
at most one, as `std::unique_ptr<ForestCombiner<L>> combiner_`, built only by
the BCF constructor. The landed virtual surface:

- `formForestResponse(f, forests, y, w)` - forest f's effective (response,
  weights) pair against the shared residual: the per-forest view a combiner
  forms so forest f's own leaf draws see a location and precision that
  reproduce the residual net of every other forest's scaled contribution. Pure
  virtual; every combiner must define it.
- `combinedFits(forests)` - the combined per-observation location over every
  forest, the pointer response_->refreshLatents/drawSigma and storeSample's
  trainingFits read. Pure virtual.
- `drawGlue(rng, sigma, y, w, forests)` / `afterCombine(forests, record,
  sampleNum, rng)` - the coupling draw and its likelihood-invariant
  post-combine move, fired at the fixed sweep points in the fixed order (a,
  aVariance, b0, b1, then the ridge rescale v). Both are inert by default (a
  combiner that only forms an additive combination need not override them).
  afterCombine's return is a REPORTING channel, not a record of whether it
  moved: each override states its own convention, 1.0 does NOT mean the state
  is unchanged, and no caller may read it that way (CORRECTED here; the base
  Doxygen is authoritative, combiner.hpp:562-572). The per-forest amplitude
  rescale returns the scale applied to the forest it reports, 1.0 if that one
  held while another travelled; the multinomial level shift returns 1.0
  unconditionally, HAVING moved, an additive move having no scale to report.
  The sweep discards the value, but `Chain::interweaveGlueRidge` - a public
  forwarder kept for the component tests - passes it through, which is how
  tests/cpp pins the ridge move's magnitude without reaching into the
  combiner's private state.
- `drawForestGlue(f, rng, forests)` - a per-forest pre-update hook, fired inside
  the sweep just before forest f's tree update with the partially updated
  forests (0..f-1 new this sweep, f..K-1 old). A no-op consuming no rng by
  default, so the additive combiners (BCF included) stay bitwise unchanged; an
  interleaved coupling (multinomial's one-vs-rest Polya-Gamma,
  docs/design/multinomial.md) draws forest f's latents here against the current
  margins, immediately before `formForestResponse(f)` reads them.
- `reportedForest()` / `testFitsAreDefined()` / `logLikelihoodIsDefined()` -
  the reporting map: which forest the scalar channels (variable counts, k,
  split probabilities) address, and whether the test-fit and log-likelihood
  channels can be computed at all. Defaults report forest 0 with both
  channels defined, so a future additive-only combiner inherits sane
  reporting for free.
- `serializeGlue(ChainStateData&)` / `restoreGlue(const ChainStateData&)` -
  glue (de)serialization into the state wire format. The base no-ops leave
  `ChainStateData::hasBCF` untouched; Chain sets it false before calling
  serializeGlue, so a combiner that owns no such flag never marks it true.
  BCFForestCombiner's serializeGlue is the sole writer of hasBCF = true - the
  serialize side owns the marker, matching the getState/setState/installForest
  contract (three glue sites, all delegating through this pair - see below).
- `setTreatment(const double*)` / `bcfGlue(double&, double&, double&) const` -
  the BCF-shaped query/mutation surface `Chain::setTreatment` and
  `Chain::bcfGlue` forward to. Both are inert/false by default; bcfGlue's
  false return is the "no glue" answer a future non-BCF combiner gives for
  free, the same answer Chain already gives when combiner_ is null.
  (BOTH RETIRED SINCE, at M4.3: `setTreatment` is gone from all four layers in
  favour of `setForestBasis(f, values, numColumns)`, and `bcfGlue` gave way to
  the ragged `totalAmplitudes`/`numForestAmplitudes`/`amplitudes` trio, of
  which bcf's three scalars are a named non-virtual reading. See
  docs/design/multiplier-combiner.md, "One mutation route".)

BCFForestCombiner<L> is static_assert-gated to `!L::hasVectorParams &&
!L::hasFunctionParams` - constant leaf only, mirroring the constraint the
two-forest Chain constructor already carried before this refactor (BCF is a
constant-leaf model end to end).

## The null-short-circuit fast path

combiner_ is nullptr for every single-forest chain, and stays the ONLY test at
every touchpoint - never a NullCombiner sentinel object, which would force a
per-sweep virtual call the single-forest chain must not pay. The landed
touchpoints (chain.hpp): setTreatment, bcfGlue, interweaveGlueRidge,
formForestResponse inside both sweep loops (run() and growForestFromRoot()),
drawGlue+afterCombine at the sweep's glue point, combinedFits() (returns the
bare `forests_[0].totalFits.data()` pointer off BCF, with no virtual call and
no copy), serializeGlue in getState, restoreGlue in both setState and
installForest, and storeSample's forest selection plus its trainingFits,
testFits, and logLikelihood branches. Every one of these is a single
`if (combiner_)` or a ternary against it - the identical shape `if (bcf_)` had
before the relocation, so the single-forest sweep is byte-identical and pays
no new indirection (equivalence 22/22 is the gate that proves it, every
commit).

Combiner methods are per-SWEEP granularity, never per-observation - the design's
"never per-observation virtual dispatch" rule (core-generalization.md) holds
here exactly as it does for ResponseModel: a combiner's O(n) loops (the
residual divide in formForestResponse, the blend in combinedFits, the glue's
Gaussian full conditionals in drawGlue) run inside one virtual call per sweep,
never inside the per-observation partition/suffstat kernels.

## BCFForestCombiner<L>: the first instance

**GENERALIZED SINCE (M4.0-M4.3, 2026-08-13 to 2026-08-14). The math this
section carried is now docs/design/multiplier-combiner.md's, and this section
is a pointer plus the combiner-hierarchy facts.** The class is no longer a
two-forest object: it is the general K-forest basis/amplitude family, each
forest contributing `m_{f,i} f_f(x_i)` with `m_{f,i} = dot(a_f, B_f(i, .))` a
contraction of that forest's own n x q_f row-major basis with its own amplitude
vector, of which bcf's `a mu + b_z tau` is the K = 2 instance
(combiner.hpp:685, chain.hpp:687-689, facade.hpp:774-776). The SPELLING is
still bcf's at every layer, which is recorded as debt in that design note.

What is combiner-hierarchy content, and stays here:

- It holds `BCFState` (the per-forest bases and their canonical flags, the flat
  ragged amplitude vector and its offsets, the per-forest amplitude priors, and
  the per-sweep combined/forestResponse/forestWeights scratch) plus a
  `const ColumnStore&` for the observation count and, in afterCombine, the data
  the ridge move touches. Built from `(data, spec, numForests)` by Chain's
  K-forest constructor.
- It is `static_assert`-gated to a constant leaf - a hierarchy constraint, not a
  model choice. The chain it is built by was Gaussian-only; M4.4 lifted that to
  gaussian, probit and logistic, the three families the calibration map has a
  latent scale anchor to state a node scale against.
- Every virtual it overrides, it overrides as ONE instance would: which
  channels it defines (`testFitsAreDefined` and `logLikelihoodIsDefined` both
  false, so storeSample NaN-flags rather than silently reporting the bare
  reported forest), which conduits it opens (`supportsResponseMutation` and
  `supportsForestWeights` both true, each with its own recorded reason), and
  which it declines (`combinedTestFits`, the counts/offset trio,
  `setActiveRows`, all left at the base's refusing default).
- `serializeGlue` is the sole writer of `hasBCF = true`, per the marker
  contract above; `restoreGlue` is a no-op on a state carrying no glue, and
  `glueIsValid` - a THIRD wire virtual, added after this refactor - refuses a
  state whose per-forest widths differ from the live ones even at an equal
  total.

Everything else - the model, the amplitude layout, the reparameterization and
its `0x1p-26` snap, the q-variate conditional and its LDL' factorization, the
general per-forest ASIS ridge (which is no longer a PROGNOSTIC-only move: its
q = 1 instance is bcf's a-move bitwise, its q = 2 fixed-variance instance the
b-move, and both are one mechanism at exponent `p = (L - q)/2`), the canonical
draw-path predicate, the single mutation route, the ragged persistence layout
and the three accumulation contracts - is docs/design/multiplier-combiner.md.

## What still re-carves when a second combiner lands

The title's premise is spent: three combining models have landed since it was
written (multinomial 2026-07-15, heteroscedastic's weight-channel variance
forest 2026-07-20, the general multiplier family 2026-08-14). It is kept as the
running record of what has and has not generalized.

The combiner's input side already generalizes beyond two forests without
having built past BCF: `combinedFits`/`formForestResponse` take the whole
forest vector, not a hardcoded mu/tau pair; `formForestResponse` already
returns a (response, weights) pair per forest, so a future variance forest can
route into the WEIGHT channel rather than the additive location - the seam
exists, unused; `drawGlue`/`afterCombine` are virtual, so each future model
owns its coupling and its internal draw order freely; `storeSample` already
asks the combiner which forest each reported channel addresses via
reportedForest() rather than hardcoding forest 0.

What did NOT yet generalize was Chain-level, not combiner-API, and was the
honest remaining work the second combiner would meet. Three of the four are now
resolved - multinomial took the first, heteroscedastic the second, the
multiplier family the fourth - leaving hurdle's the only one standing:

- The combined-fit OUTPUT was a single n-vector: `combinedFits` returned
  `const double*` and `results.trainingFits` was one channel. RESOLVED by
  multinomial: its n x K combined object was carried by a location-count seam
  (`numReportedLocations()`, docs/design/multinomial.md) widening combinedFits
  and the training/test writes; `refreshLatents`/`drawSigma` did NOT widen
  (softmax needs no per-observation location on the response side, its K PG
  draws living in the combiner's drawForestGlue), so the predicted
  three-call-site change was really one seam.
- Chain holds exactly one `response_` and one `sigma_`. Heteroscedastic's
  variance forest needed either the unused weight-channel route above or a
  per-observation sigma - a decision this plan explicitly deferred. RESOLVED,
  2026-07-20 (docs/design/heteroscedastic.md): the WEIGHT-channel route was the
  right one, and it landed as a Chain-side nullable `varianceForest_`
  (chain.hpp:680-682) rather than through a combiner at all. The Chain
  constraint itself stands unchanged; what is settled is the decision, not the
  constraint. The combiner's own weight-channel seam is therefore STILL the
  unused route a future combiner-HOSTED variance forest would take.
- Chain holds a single-leaf-type `std::vector<Forest<L>>` (one L for every
  forest). Hurdle's per-forest response families (an occupancy forest under
  one family, a positive-part forest under another) break that invariant;
  it is not a property the combiner API constrains. STILL OPEN, and see
  "Anticipated" below for why hurdle no longer asks for it.
- `ChainStateData`'s glue fields were BCF-shaped (a/aVariance/b0/b1). RESOLVED
  by the multiplier family (M4.3): they are now `hasBCF`, `amplitudeWidths`,
  `amplitudes` and `amplitudeVariances` - RAGGED, with the widths travelling
  because a TOTAL IS NOT A LAYOUT, `q = (1, 3)` and `q = (2, 2)` both carrying
  four amplitudes (combiner.hpp:92-102) - and the four named scalars survive
  only as a hand-written K = 2 reading, non-authoritative (:103-109). The
  bullet's CONCLUSION is what the arc vindicated and it stands: a non-BCF
  combiner overrides `serializeGlue`/`restoreGlue` rather than reaching for an
  accessor, so the interface point was already right. It need not always write
  anything: the multinomial combiner serializes NOTHING (docs/design/
  multinomial.md), redrawing its per-sweep Polya-Gamma latents against the
  restored forests structurally, so restore is structural, not bitwise. The
  interface point did grow a THIRD virtual, `glueIsValid` (combiner.hpp:675,
  :1018-1028), the layout check `stateIsValid` routes through so a same-total
  different-layout state cannot be written through the live offsets.

### What the multiplier family closed, and what it opened

The combiner API's own input side is now general in K AND in per-forest basis
width: `setForestBasis` (combiner.hpp:511-517), the amplitude trio
`totalAmplitudes`/`numForestAmplitudes`/`amplitudes` (:550-560), and a ragged
glue wire block (:666-675). The math is docs/design/multiplier-combiner.md's
and is not restated here.

**What the virtual surface grew, and the general rule it grew under.** Since
2026-07-14 the base gained `combinedTestFits`, `setForestBasis`,
`supportsCountsMutation`/`setCounts`/`setCategoryOffset`/
`setCategoryTestOffset`, the amplitude trio, `forestReportingIsDefined`,
`numReportedLocations`, `numVariableCountForests`/`variableCountForest`,
`supportsResponseMutation`, `supportsForestWeights`, `setActiveRows` and
`glueIsValid` (combiner.hpp:487-676). Every one of them follows one rule, and
it is the sentence worth landing: **each is a CAPABILITY predicate defaulting
to the REFUSING answer, so a future combiner stays refused at the bridge until
it is audited** (combiner.hpp:526-528, :634-638, :645-647). And never a
forest-count test, because a K-forest multinomial defeats one - which is why
the bridge probes `totalAmplitudes() != 0` instead (C_interface.cpp:764-770).

What still does NOT generalize, after M4:

- (i) The single-leaf-type forest vector, unchanged (third bullet above).
- (ii) Non-Gaussian. The K-forest constructor hardcodes `GaussianResponse` and
  `family_ = ResponseFamily::gaussian` (chain.hpp:702-705), and
  `createBCFSampler` carries a single `SamplerFacade<ConstantGaussianLeaf>`
  instantiation (facade.hpp:786-788). That is M4.4.
- (iii) The NAMING. The general family is spelled `BCFSpec` /
  `BCFForestCombiner` / `ChainStateData::hasBCF` / the `"bcf"` state block, so
  every layer reads "BCF" where it means "carries amplitudes". Recorded as debt
  (root TODO `bcf-naming-generalization`, described in
  docs/design/multiplier-combiner.md); its mitigation is the capability-probe
  rule above, which keeps no consumer-visible answer keyed on the name.

## Anticipated (multinomial now built)

- Multinomial: now BUILT (docs/design/multinomial.md), the second combiner and
  the extraction's occasion. K SYMMETRIC forests coupled through a softmax
  likelihood with an INTERLEAVED one-vs-rest Polya-Gamma augmentation:
  drawForestGlue(k) draws omega_k against the current margins just before forest
  k's update, formForestResponse(k) forms category k's PG working response
  against its log-sum-exp margin, and afterCombine runs a level-centering move
  (not a category-conditional Gaussian glue). It took the combined-fit n x K
  widening named above (the numReportedLocations() location-count seam) but no
  wire-format bump - the softmax combiner serializes nothing.
- Heteroscedastic (HBART): a mean forest plus a variance forest with
  multiplicative-positive leaves. The variance forest's natural route is the
  WEIGHT channel formForestResponse already returns, unused; a non-integrable
  leaf and its own MoveStrategy are out of this refactor's scope entirely, as
  is the per-observation-sigma-vs-weight-channel decision above.
  (SUPERSEDED 2026-07-19, docs/design/heteroscedastic.md: the variance leaf is
  CONJUGATE, not non-integrable - a scaled-inverse-chi-squared scale leaf reusing
  the existing conjugate move, no new MoveStrategy. The WEIGHT-channel route was
  correct. See that note for the full design.)
- Hurdle: a binary-occupancy forest plus a positive-part forest sharing
  predictors through the data handle. Response families differing per forest
  breaks Chain's single response_ - Chain-level, not combiner-API, per above.
  (LANDED 2026-07-20 R-SIDE, docs/design/hurdle.md: the two parts are
  conditionally independent, so hurdle composes two ordinary fits glued at
  report time and this Chain break stays UNBUILT; the engine two-response route
  awaits a genuinely coupled model - zero-inflation or Heckman selection.)
- Grouped-x-multi-forest: GroupedResponse decorates the response chain BELOW
  combining - the sweep feeds the combiner `y = response_->workingResponse()`,
  already group-adjusted if response_ is a GroupedResponse - so a grouped
  decorator and a combiner compose without either knowing the other. This is
  the composition architecture-numerical-review.md's debt #1 flagged as
  impossible while BCF was a hardcoded Chain special case; this refactor makes
  it expressible. It does not build grouped-BCF or any grouped-multi-forest
  model.

## Post-mutation fit rebuild was forest-0-only (CLOSED)

`revalidateTrees`, `rebuildFitsFromParameters`, `applyNewData`, and
`recoverTreeParameters` (chain.hpp) all read `forests_[0]` directly and did not
loop over forests_ the way the sweep, storeSample, getState/setState, and
installForest already did. This was a pre-existing forest-split interim
(forest-split-bcf.md) this refactor did not close: a BCF chain's tau forest
could not ride setPredictor's whole-data-replacement transaction or a
warm-start's tree revalidation. Nothing about the combiner blocked widening
them to loop over forests_; it was simply out of this plan's scope.

CLOSED, 2026-08-10 to 2026-08-12, by the multiforest-predictor-mutation arc
(S1-S4; recorded at docs/design/model-space-survey.md:383-400). Every one of
those paths now loops the forests, including the heteroscedastic variance
forest, and the acceptance rule resolved to PER-SAMPLER: a row installs only if
it empties no leaf in any tree of any forest of any chain.
`refuseMultiForestTransactionalUpdate` is retired from all four transactional
entries, and feature-matrix.md carries `bcf x setPredictor = S`.

## A standalone hierarchy, not a ResponseModel subclass

ForestCombiner<L> is its own hierarchy beside ResponseModel, not a subclass or
member of it, even though architecture-numerical-review.md's debt #1 named the
gap as combining belonging "on the ResponseModel side." VD confirmed reading
that as a responsibility statement, not a class-hierarchy instruction
(docs/plans/forest-combiner.md, Open questions): combining is a response-side
CONCERN, owned by a dedicated object, not necessarily a ResponseModel
subclass. The one-line why: ResponseModel's interface is per-observation-
location (workingResponse/workingWeights, refreshLatents, drawSigma - one
scalar per observation, read by kernels that never see a forest); a
combiner's job is per-forest residual formation over the WHOLE forest vector,
a materially different shape. Folding combining into ResponseModel would force
forest-vector methods and the leaf-type template onto every single-forest
response family (Gaussian, probit, logistic, AFT) to serve a case none of
them need. Keeping the two orthogonal - Chain feeds the combiner
`y = response_->workingResponse()`, already whatever response_ decorates it
into - is exactly what lets a combiner and a grouped decorator compose without
a class per pairing, the point the "anticipated, not built" section above
relies on.

## Header location

Defined in chain.hpp beside Forest<L>, not a separate combiner.hpp: the ridge
move (afterCombine) reaches Forest<L>'s tree/fit buffers and saved-tree
FlatNodes directly, so a combiner.hpp would either fight chain.hpp's include
order or force a Forest<L> extraction this neutral refactor had no reason to
carry. Extraction to src/bartcore/combiner.hpp DID land with the second combiner
(multinomial, docs/design/multinomial.md): Forest<L>, ForestResponse,
ForestCombiner<L>, BCFForestCombiner<L>, and the serializable state/spec structs
moved there as pure motion and chain.hpp includes it - the second consumer
shaped the split, as intended, rather than a guess made against one case.

## Verification

Landed under two bitwise gates, run at every commit. The existing
single-forest anchor (benchmarks/R/equivalence.R vs
equivalence-ac6ec2c.rds, 22/22 IDENTICAL) guards the fast path; it says
nothing about BCF, since equivalence fits only through dbarts(), which never
builds a multi-forest sampler. benchmarks/R/bcf-equivalence.R is the BCF
analogue this plan added (step 1): five scenarios (default two-forest, a
restricted moderator, the updateA/updateB toggle, a WEIGHTED fit, and a
setTreatment mutation), each recording six channels (both forests' raw fits,
the glue, sigma, and the reported train/varcount channels) and asserting
bitwise identity, not tolerance. bcf-exact{,-restricted,-weak}.R (a
closed-form posterior match to Monte Carlo error) backstops correctness the
bitwise fixture can't speak to - a wrong-but-neutral relocation would pass the
fixture and fail there. tests/cpp gained a BCF growForestFromRoot case (that
sweep branch is unreachable from R) and tightened testBCFInterweaveKeepTrees
to numThin > 1 (the saved-slot addressing the ridge rescale touches was
previously pinned only at numThin = 1), closing the two component-gate holes
the R-level fixture cannot reach on its own.

## Landing notes

Commit 1 = 7145e29. Record a bitwise BCF fixture and close the component gate
holes (test-only, no engine change). benchmarks/R/bcf-equivalence.R added in
equivalence.R's record/compare idiom: five scenarios (default, restricted
moderator, updateA/updateB toggle, weighted, setTreatment) x six channels
(bartcoreForestFits for both forests, bartcoreBCFGlue, sigma, result$train,
result$varcount), asserted with identical(), plus a settings guard that
refuses a compare against a baseline recorded under different knobs.
benchmarks/baselines/bcf-equivalence-99205ee.rds recorded from the
pre-refactor HEAD. tests/cpp gained testBCFGrowForestFromRoot (the
growForestFromRoot sweep's BCF branch, unreachable from R) and tightened
testBCFInterweaveKeepTrees with a numThin = 3 block pinning the thinned
saved-slot addressing the ridge rescale touches. Gate: equivalence 22/22 +
full tinytest unchanged + tests/cpp (grown) + the new fixture identical to
itself.

Commit 2 = 691697b. Introduce the forest combiner and relocate the combining
math (engine, byte-identical). ForestCombiner<L> and BCFForestCombiner<L>
introduced; combiner_ replaces bcf_; forestMultiplier/formForestResponse/
combinedFits moved into BCFForestCombiner. DEVIATION recorded: renaming the
owner forced every glue reader (storeSample, getState/setState/installForest,
setTreatment, bcfGlue) to be mechanically rewritten to read through combiner_
in this step rather than step 4, since there was no other byte-identical way
to keep the no-mu/tau-hardcoding contract while bcf_ no longer existed;
storeSample's combinedFits() call was pulled forward here for the same
reason. The rewrite went through a transitional `glueState()` virtual
accessor (returning `BCFState*`/`const BCFState*`) added to ForestCombiner<L>
so Chain's still-resident drawGlue/interweaveGlueRidge/storeSample/getState/
setState bodies could reach BCFForestCombiner's internal glue struct directly;
this was acknowledged as transitional and fully retired in commit 4, once
drawGlue/afterCombine/serializeGlue/restoreGlue took over its call sites.
Gate: equivalence 22/22 + BCF fixture identical + tests/cpp + bcf-exact.R
quick + full tinytest.

Commit 3 = b8b0b4c. Move the coupling draw and ridge move into the combiner
(engine, byte-identical). drawGlue and interweaveGlueRidge (the ridge rescale)
became virtual ForestCombiner methods, called from Chain at the same sweep
points and in the same rng order as before. DEVIATION recorded: afterCombine
returns double (then the applied ridge scale c, 1.0 when the move was skipped
or inert - see the corrected LIVE contract at "The ForestCombiner<L>
hierarchy" above, which this historical sentence predates: 1.0 no longer means
unchanged) specifically so the cpp component tests can read it through
Chain::interweaveGlueRidge, the public forwarder kept for exactly this; the
sweep itself discards the return value. Doxygen comments relocated with the
math were redirected from docs/plans/bcf-ridge-interweaving.md (and a
"docs/plans Status" pointer) to docs/design/bcf.md, since that plan's
landing had already folded into bcf.md's Burn-in and Landing sections - a
relocated comment should point at where the math is documented going
forward, not at the one-time plan that shipped it. Gate: equivalence 22/22 +
BCF fixture identical + tests/cpp + bcf-exact.R quick+full +
bcf-exact-restricted.R + bcf-exact-weak.R + tinytest.

Commit 4 = c90ecf2. Route reporting, serialization, and treatment through the
combiner (engine, byte-identical). storeSample's forest selection and channel
definedness, getState/setState/installForest's glue (de)serialization, and
setTreatment/bcfGlue all moved onto the virtual surface
(reportedForest/testFitsAreDefined/logLikelihoodIsDefined/serializeGlue/
restoreGlue/setTreatment/bcfGlue); ChainStateData.hasBCF is now set by
Chain before calling serializeGlue and owned on the serialize side.
No deviation from the plan text: the transitional glueState() accessor from
commit 2 was retired in full (every call site now goes through the named
virtuals), and Chain no longer names BCFState anywhere. Gate: equivalence
22/22 + BCF fixture identical + tests/cpp (incl. the state fuzzer) +
bcf-exact.R quick + full tinytest 2832, no regen.

Commit 5 = this commit (docs). docs/design/forest-combiner.md added (this
file); docs/plans/multi-forest-models.md's forest-combiner blocker marked
discharged; docs/plans/architecture-numerical-review.md's debt #1 marked
closed; the repo-root TODO's forest-combiner entry retired and the
multi-forest-models entry's prerequisite clause updated. Gate: R CMD check
man unaffected (no man/ touched); full tinytest 2832 PASS, 0 FAIL (docs-only,
no growth from commit 4).

Gates held at every commit: equivalence 22/22 IDENTICAL vs
equivalence-ac6ec2c.rds throughout (the single-forest fast path never moved);
the BCF fixture's five scenarios x six channels identical at every commit
from its commit-1 recording forward; bcf-exact.R quick and full, plus
bcf-exact-restricted.R and bcf-exact-weak.R, unmoved at every BCF-touching
commit - the full run's 2b E[mu] gap held at 0.0172, matching the
pre-refactor record (forest-split-bcf.md); full tinytest held at 2832 with no
snapshot regeneration across all five commits; the tests/cpp component suite,
grown at commit 1, passed at every subsequent commit.

Confirmatory bench (quiet window, arm64, single-forest hot path).
bench-sampler compare vs benchmarks/baselines/bench-sampler-4008675.csv at
the docs commit: "OK: no metric regressed more than 5%", zero flags - run
arms 1.018-1.043x baseline (within the recorded sub-5% noise band),
setPredictor arms 0.987-1.012x. The null-short-circuit fast path costs the
single-forest sweep nothing, as designed. The plan is CLOSED.
