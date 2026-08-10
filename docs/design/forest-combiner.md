# Forest combiner: design

Status: LANDED, 2026-07-14. Promotes the BCF glue that lived as a Chain-side
special case (bcf_, drawGlue, combinedFits, formForestResponse,
forestMultiplier, the if(bcf_) sweep branches - the shape bcf.md's Forest
split shipped in) into a polymorphic ForestCombiner<L> hierarchy beside
Forest<L> (src/bartcore/combiner.hpp since the multinomial extraction below;
src/bartcore/chain.hpp at this refactor's own landing). BCFForestCombiner<L> is
its first instance, the math it carries unchanged from bcf.md and
bcf-ridge-interweaving's landing; MultinomialForestCombiner<L>
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
  afterCombine returns the scale its move applied (1.0 when it makes none);
  the sweep discards the return value, but `Chain::interweaveGlueRidge` - a
  public forwarder kept for the component tests - passes it through, which is
  how tests/cpp pins the ridge move's magnitude without reaching into
  BCFForestCombiner's private state.
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

Holds `BCFState` (the borrowed treatment vector z, the glue a/b0/b1/aVariance,
the priors, the updateA/updateB switches, and the per-sweep combined/
forestResponse/forestWeights scratch) plus a `const ColumnStore&` for the
observation count and, in afterCombine, the data the ridge move touches. Built
from `(data, spec)` by Chain's BCF constructor. The math is unchanged from
bcf.md and bcf-ridge-interweaving's landing, relocated verbatim:

- `formForestResponse` divides the residual by forest f's scale multiplier (a
  for mu, b_z for tau) and scales the weight by the multiplier squared - the
  resid/m, w*m^2 pair that reproduces bcf.md's model equation (y = a mu +
  b_z tau + eps) in each forest's own constant-leaf node sums. A multiplier
  indistinguishable from zero at `sqrt(DBL_EPSILON)` (2^-26) snaps both to
  exactly 0.0 instead of dividing by a floored near-zero value (2026-08-10,
  docs/plans/zero-weight-exactness.md); it composes with an optional
  caller-settable per-forest weight `s_{f,i}` installed by
  `Chain::setForestWeights`, applied as one further multiplicative factor on
  the weight channel right after this call returns, before the tree loop.
- `combinedFits` blends `a * mu + b_z * tau` per observation.
- `drawGlue` draws the Gaussian full conditionals for a (with its half-Cauchy
  scale-mixture auxiliary aVariance) and b0/b1 (docs/design/bcf.md).
- `afterCombine` is the interweaving (ASIS) rescale of the prognostic glue
  ridge (docs/design/bcf.md, "Burn-in" and bcf-ridge-interweaving's folded-in
  history), reaching Forest<L>'s treeFits/totalFits/test fits and saved-tree
  FlatNodes to apply the same scale c the glue a absorbed.
- `reportedForest()` returns 0 (the prognostic forest); `testFitsAreDefined()`
  and `logLikelihoodIsDefined()` both return false, since BCF's API carries no
  test treatment vector to blend an off-sample a*mu + b_z*tau, and the blended
  per-observation location is not visible to the response model to score -
  storeSample NaN-flags both channels rather than silently misreporting the
  bare prognostic forest.
- `serializeGlue` sets hasBCF = true and copies a/aVariance/b0/b1 into
  ChainStateData; `restoreGlue` is a no-op when the state carries none, so a
  mismatched restore leaves the glue at its constructed values rather than
  clobbering them with zeros.

## What still re-carves when a second combiner lands

The combiner's input side already generalizes beyond two forests without
having built past BCF: `combinedFits`/`formForestResponse` take the whole
forest vector, not a hardcoded mu/tau pair; `formForestResponse` already
returns a (response, weights) pair per forest, so a future variance forest can
route into the WEIGHT channel rather than the additive location - the seam
exists, unused; `drawGlue`/`afterCombine` are virtual, so each future model
owns its coupling and its internal draw order freely; `storeSample` already
asks the combiner which forest each reported channel addresses via
reportedForest() rather than hardcoding forest 0.

What does NOT yet generalize is Chain-level, not combiner-API, and was the
honest remaining work the second combiner would meet; multinomial (below) has
since resolved the first and last of these, leaving heteroscedastic and hurdle
the middle two (recorded here so multi-forest-models.md plans against reality,
not against what the API merely gestures at):

- The combined-fit OUTPUT was a single n-vector: `combinedFits` returned
  `const double*` and `results.trainingFits` was one channel. Multinomial's
  n x K combined object was carried by a location-count seam
  (`numReportedLocations()`, docs/design/multinomial.md) widening combinedFits
  and the training/test writes; `refreshLatents`/`drawSigma` did NOT widen
  (softmax needs no per-observation location on the response side, its K PG
  draws living in the combiner's drawForestGlue), so the predicted
  three-call-site change was really one seam.
- Chain holds exactly one `response_` and one `sigma_`. Heteroscedastic's
  variance forest needs either the unused weight-channel route above or a
  per-observation sigma - a decision this plan explicitly deferred (TODO,
  architecture-numerical-review.md).
- Chain holds a single-leaf-type `std::vector<Forest<L>>` (one L for every
  forest). Hurdle's per-forest response families (an occupancy forest under
  one family, a positive-part forest under another) break that invariant;
  it is not a property the combiner API constrains.
- `ChainStateData`'s glue fields are BCF-shaped (a/aVariance/b0/b1). A non-BCF
  combiner overrides `serializeGlue`/`restoreGlue` rather than reaching for an
  accessor - the interface point is already right. It need not always write
  anything: the multinomial combiner serializes NOTHING (docs/design/
  multinomial.md), redrawing its per-sweep Polya-Gamma latents against the
  restored forests structurally, so restore is structural, not bitwise. Only a
  combiner carrying un-recoverable scalar glue of its own shape would need a
  wire-format bump (the flat-format-v2 scheme public-surface.md anticipates)
  first - not every multi-forest model does.

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

## Post-mutation fit rebuild stays forest-0-only

`revalidateTrees`, `rebuildFitsFromParameters`, `applyNewData`, and
`recoverTreeParameters` (chain.hpp) all read `forests_[0]` directly and do not
loop over forests_ the way the sweep, storeSample, getState/setState, and
installForest already do. This is a pre-existing forest-split interim
(forest-split-bcf.md) this refactor did not close: a BCF chain's tau forest
cannot yet ride setPredictor's whole-data-replacement transaction or a
warm-start's tree revalidation - both are single-forest paths today, silently
scoped to mu. Nothing about the combiner blocks widening them to loop over
forests_; it was simply out of this plan's scope (it touches no combining
math and no per-sweep behavior).

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
returns double (the applied ridge scale c, 1.0 when the move is skipped or
inert) specifically so the cpp component tests can read it through
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
