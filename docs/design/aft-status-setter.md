# An AFT censoring-status setter, and the SBC arms it enables

Status: LANDED - slice 1 (section 8), 2026-09-07 (fcd60feb, e20c6462, f9bc9260), slice 2, 2026-09-07 (5125f7cd,
2c766437), slice 3's prior-draw entry and variance-surface accessor, 2026-09-13 (31149f5b), and slices 3 and 4's
heteroscedastic SBC arms, 2026-09-13 (4424e69a). Every slice has landed, and all three arms - aft and the two
heteroscedastic ones - are in the SBC matrix.

Amended by [pure-c-header](../plans/pure-c-header.md#pure-c-header): the flat C header creates no sampler and
no longer declares the predictor, test-data, weight, active-row, per-forest, state,
tree-extraction or augmentation entries - each is a method on the R sampler object the
handle is now read from. The `retired:` cites below name constructs that are gone; what
this record says about the R and engine sides still holds.

Lets a live `family = "aft"` sampler take a new per-observation censoring status, so the censoring structure stops being
fixed at creation. The enabled item is SBC coverage: aft is out of the matrix for exactly this reason
([Decision - scope](../plans/sbc-family-tiers.md#decision---scope)), and the heteroscedastic composition inherits the gap
([6. Gates](aft-variance-forest.md#6-gates)).

## 1. What is fixed today

[`AFTResponse`](../../src/bartcore/model.hpp) derives `censoredIndices_` and `censorBound_` once, in its constructor, from
the status the bridge reads off the control attribute ([`applySurvivalAttribute`](../../src/R_interface_bartcore.cpp)),
keeping no pointer to it. [`AFTResponse::setResponse`](../../src/bartcore/model.hpp) replaces the log-times and refreshes
each bound BY INDEX, so a new response reaches the old censoring structure. Whole-data replacement is refused outright
(["fix the censoring structure at creation"](../../src/R_interface_bartcore.cpp)).

Two facts make a setter cheap. `censorBound_` shadows the OBSERVED log-time at every censored row, so that data survives the
latent overwrite: row `i`'s observed time is `censorBound_[k]` where censored and `logT_[i]` where not. And the contained
[`GaussianResponse`](../../src/bartcore/model.hpp) takes its posterior degrees of freedom from the positive-weight count,
not the event count ([`GaussianResponse::sigmaDegreesOfFreedomForTesting`](../../src/bartcore/model.hpp)), so a status
change moves the sigma posterior only through the imputed values. Creation validates three things - a real vector, length n,
elements 0 or 1 (NaN fails both) - with no minimum event count; an all-censored status is accepted, sigma then identified
through the truncation and the prior.

## 2. Surface

**A, a dedicated `$setSurvivalStatus(status)`.** One method, one bridge entry, no other family's signature touched. But it
splits an operation the engine does as one, and alone it carries section 3's rebuild trap.

**B, widen `$setResponse(y, status = NULL)`.** `status` goes LAST, name-only and no prefix of `updateScale` or
`updateState`, so no positional call changes meaning and the warning
(["setResponsePositionalUpdateScale"](../../R/dbarts.R)) cannot fire on it. `NULL` is today's behaviour exactly; non-aft
families refuse a non-null `status` by name. A status-only change is `s$setResponse(s$data@y, status = s2)`, `data@y` being
the log time [`bartcoreSamplerSetResponse`](../../R/bartcore.R) mirrors every later write into.

**C, both.** Two entries for one engine operation, two refusal sets, two Rd blocks. RECOMMEND B: the transaction the engine
wants and the one the SBC harness issues, with A's door open at the cost of one forwarding method. Semantics, inherited
rather than invented:

- **Transaction.** Type, length and value checks run in the bridge before the engine is touched; the R5 method mirrors into
  `attr(control, "bartcore.survival")` only after the `.Call` returns. Validate before mirror; the R5 `$setWeights`
  `tryCatch` rollback is unnecessary, a refusal installing nothing. The R5 method coerces with `as.double`, leaving the type
  refusal to the handle.
- **updateScale.** The status re-anchors nothing, and a heteroscedastic aft takes either flavor - a re-anchoring swap
  restates the variance forest's scale leaf and surface ([`Chain::reanchorVarianceForest`](../../src/bartcore/chain.hpp)) -
  with both R entry points defaulting it FALSE. Not an invariant: the transform is over `logT_`, which equals the observed times only just after a `setResponse` -
  [`AFTResponse::setOffset`](../../src/bartcore/model.hpp) at `updateScale = TRUE` re-anchors on latents.
- **Latents.** The redraw stays, matching [`AFTResponse::setResponse`](../../src/bartcore/model.hpp) and the probit pattern:
  event to censored takes the row's own observed time as its new bound and redraws above it; censored to event restores that
  time, which is data, not a draw. A caller driving one sweep per call needs the imputed value, the next sweep's mean forest
  reading the working response before [`AFTResponse::refreshLatents`](../../src/bartcore/model.hpp).
- **Caveat, masked samplers.** The memcpy is mask-blind, the redraw is not, so the status-only idiom resets EVERY inactive
  censored row's latent to its bound and leaves it there until the mask clears. Pre-existing `setResponse` behaviour, newly
  advertised; documented on the method.

## 3. Engine

The interface is a fork. **(i) Widen the pure virtual:** `ResponseModel::setResponse` is `= 0` with eight concrete overrides
(gaussian, probit, ordinal, logistic, multinomial, aft, t, nbinom), and a default argument does not help, being statically
bound; all eight change signature, as do `Chain`, `Sampler`, [`SamplerBase`](../../src/bartcore/facade.hpp) and
`dbarts_sampler_setResponse`'s site in src/C_interface.cpp. **(ii) A separate non-pure `ResponseModel::setSurvivalStatus`,**
default no-op with an `AFTResponse` override - the [`ResponseModel::setVarianceSurface`](../../src/bartcore/model.hpp)
shape - called by the widened `bartcore_setResponse` BEFORE `setResponse`. RECOMMEND (ii): it edits one response model
instead of eight, and ordering it first gives the joint call exactly one latent redraw and one working-response rebuild; the
`setResponse` facade virtual and src/C_interface.cpp stay untouched (landed: the state handshake below needs its own digest
and reapply pair, so `AFTResponse` ends up with three overrides - `setSurvivalStatus`, `survivalDigest`,
`reapplySurvivalStatus` - not one, each a new [`SamplerBase`](../../src/bartcore/facade.hpp) virtual on
`weightsDigest`/`reapplyWeights`'s shape).

The override rebuilds `censoredIndices_` and `censorBound_` from the status and the observed times in force -
`censorBound_[k]` where currently censored, `logT_[i]` where not - into temporaries swapped in once the pass succeeds, so it
is correct alone; a joint call's memcpy then supersedes the bounds with the new y. Three load-bearing constraints:

- **Status first.** The reverse order builds newly censored rows' bounds out of latents - wrong data, silently. Neither call
  can fail after the bridge validates, which refuses a non-null status off aft by family as
  [`bartcore_setData`](../../src/R_interface_bartcore.cpp) does, so the pair is one transaction. `Sampler` fans the status
  to every chain, each rebuilding its own structure and redrawing off its OWN generator, the `setWeights` rule.
- **In place.** `logT_` is memcpy'd, never swapped or assigned: the contained Gaussian borrows `logT_.data()` as its own
  response. The one `assign`, [`AFTResponse::setData`](../../src/bartcore/model.hpp), re-hands the pointer; the temporaries
  above are the two index vectors only.
- **The rebuild.** `refreshLatents` early-returns at an empty censored set BEFORE `rebuildWorking`, so a bare setter ending
  in the redraw leaves `yRescaled_` holding the old censored rows' latents after a flip to all events. Under (ii) the
  `gaussian_->setResponse` that follows rebuilds it unconditionally - an argument for the joint surface, and the trap a
  future A must handle.

Three compositions need no code. The **mask**: the index rebuild and the observed-time restore are mask-blind, membership
and censoring being independent per-row facts, while the redraw inherits `refreshLatents`'s skip, section 2's caveat. The
**variance surface**: `variance_` borrows a per-row vector reallocated only from the whole-data arm aft refuses, so the
pointer stays live ([`Chain::installVarianceSurface`](../../src/bartcore/chain.hpp)). The **sigma posterior**: section 1's
degrees-of-freedom fact makes it invariant, checkable rather than argued (gate (b)).

## 4. State

Two changes, neither a new block nor a format-version move.

**Restore by index.** [`AFTResponse::restoreLatents`](../../src/bartcore/model.hpp) is today an unconditional memcpy over
all n rows. Once the status can move that loses data permanently: a row censored when the state was stored and an event
after the flip has its OBSERVED time overwritten by the donor's latent draw, unrecoverably - `refreshLatents` walks
`censoredIndices_` only, `censorBound_` no longer holds the row, and
[`AFTResponse::computeLogLikelihood`](../../src/bartcore/model.hpp) scores it at the fabricated time thereafter. The path is
the DEFAULT one: `dbartsControl` defaults `updateState = TRUE` so `$run` stores state while the mutators store only on an
explicit `TRUE`, so run, set status, save, load hands [`getPointer`](../../R/dbarts.R) a control with the new status and a
state from before it. So restore the censored indices rather than the whole vector: an event row's log-time is data a state
has no business overwriting. The RESTORE CONTRACT holds (a write plus a working rebuild reading neither sigma nor the
surface), and section 3's reconstruction becomes valid everywhere.

**A status handshake, the weights.digest shape.** `bartcore_storeState` writes a `survival.digest` attribute beside
["weights.digest"](../../src/R_interface_bartcore.cpp): top level, the status fanning to every chain; eight raw bytes;
ABSENT is not a mismatch. `bartcore_setState` compares it against the destination's live digest, which `AFTResponse` derives
from `censoredIndices_` with no status vector stored, and on a difference reapplies the structure: bounds rebuilt from the
reconstruction, the censored set redrawn off each chain's own restored generator. That is `reapplyWeights`'s
reconciliation: silent, deterministic, a measured no-op for a family carrying no status.

The mirror onto `attr(control, "bartcore.survival")` stays, being what re-creation reads: `getPointer` rebuilds from the
stored triple, so the re-created sampler takes the CURRENT status, as `data@y` and `data@counts` do. The handle has no
mirror - `bartcoreSampler` returns an environment holding a pointer and the predictors - so its wrapper ships a weaker
persistence contract and gate (d)'s check is R5-only. Residue after both, accepted: a row censored on both sides whose
restored latent sits below a bound a y change moved, healed by the next sweep's refresh before sigma is drawn
([`Chain::run`](../../src/bartcore/chain.hpp)).

## 5. Flat C API

R bridge only: the shipped header gains no entry beside
retired: [`dbarts_sampler_create`, `DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h), so the hash does not move and no
`LinkingTo` consumer recompiles. None wants it - stan4bart on bartcore maps an `"aft"` token but fits gaussian and probit
only, calling no bridge entry, and bartCause on dbarts-1.0 has no survival code. The door is one X-macro entry shaped like
retired: [`dbarts_sampler_setActiveRows`](../../inst/include/dbarts/dbarts.h), a capability-status `int` return refusing off aft;
adding it re-bakes the hash, pre-release a note rather than a migration. One consumer-visible message does change:
`bartcore_setData`'s aft refusal ("fix the censoring structure at creation") becomes false as stated. The refusal STAYS - a
whole-data replacement may change n, which the status is stated over - but the reason is restated to name the channel that
now serves the caller.

## 6. SBC

The aft arm REBUILT its fit every replication, alone among the arms, which is why it needed an anchor leaf scale named as
its own `node.prior` and an offset zeroing each rebuild's `prior.mean`: the two pinned the leaf prior and the shift against
a transform every rebuild re-derived from `range(y0)`. Both went with the conversion below. It joined
[`sbcMatrixConfigs`](../../benchmarks/R/sbc.R) and the workflow matrix
(["config: gaussian"](../../.github/workflows/sbc.yaml) was the nearest arm) only once its ladder and its R = 200 run were
read; the admission is in the landing note below.

**The aft arm.** The setter converts it to the reused-sampler shape every other arm has: one sampler built once at a fixed
build response through `dbarts(x, cbind(time, status), family = "aft")`; per replication a prior draw of `(f, sigma)`, then
`logT0 = f0 + sigma0 * eps`, `status0 = logT0 <= logC` against the censoring times
[`sbcAddCensoring`](../../benchmarks/R/sbc.R) pins, `y0 = pmin(logT0, logC)`, an overdispersed second prior draw,
`$setResponse(y0, status = status0, updateScale = FALSE)` and one `$run`. The transform never re-anchors, so both pins go.
Nine functionals: `avg.f`, `sigma`, `f.star` at the five test points, and two only this family has -
`S(t0 | x*) = 1 - Phi((log t0 - f(x*)) / sigma)`, the reported deliverable, at the first test point and `t0 = 1`, the build
response's own median survival time, which centres the ratio at 0 under the prior draw rather than in a tail S would
eventually underflow to an atom in; and `logT0[i]` at the LOWEST-INDEXED censored row, ranked against that row's posterior
latents read with `getLatents` by the per-sample `run(0, 1)` idiom the discrete functionals use. The censored set is random
per replication and can be empty; such a replication contributes no rank there and that functional's R is reported
separately. It alone ranks the truncated-normal imputation, whose channel [`sbcCheckAftLatents`](../../benchmarks/R/sbc.R)
gates beside the other arms' wiring checks: an event row's latent is its observed log time exactly, a censored row's sits
strictly above its bound.

Budget: measure a burn ladder first, which needs an aft branch in [`sbcFamilySpec`](../../benchmarks/R/sbc.R) - draw, fit,
burnRun, sample - and not only an entry in [`sbcFamilyConfig`](../../benchmarks/R/sbc.R), since
[`sbcBurnLadder`](../../benchmarks/R/sbc.R) dispatches through the former. That branch IS the reuse conversion, and it lands
with no burn pre-registered: [`sbcBurnSweeps`](../../benchmarks/R/sbc.R) carried NA for the arm, so the verdict run refused
by name until `burn-aft` was read, as the two latent BCF arms did. A sweep is a gaussian one plus a truncated-normal draw per
censored row; at a gaussian-like cost and a t-like burn, R=200, L=150, thin=30 lands near 5-10 minutes, hence a 20-30 minute
timeout at the workflow's ~3x rule. The arm's functionals raise
[`sbcMatrixFunctionals`](../../benchmarks/R/sbc.R) and widen the Bonferroni'd band for every arm, which stales the M = 30
prose in the workflow and the harness comment beside the constant; a wider band can only turn a FLAG into a PASS, so
recorded verdicts stay readable.

**The heteroscedastic arm.** None exists, the plan's reason being that the prior-draw entries are FOREST-ONLY
([`Chain::sampleTreesFromPrior`](../../src/bartcore/chain.hpp),
[`Chain::sampleNodeParametersFromPrior`](../../src/bartcore/chain.hpp) both leave the variance forest untouched, a
documented contract), so `s(x)` has no prior-draw path. Three ways to get one:

- **A, widen the two existing entries.** Contract-breaking and stream-shifting far beyond SBC: `sampleTreesFromPrior` is
  `bart2`'s default init.
- **B, draw the variance trees in R and install them through `setState`'s `variance.*` blocks.** No engine change, but it
  restates [`ConstantVarianceLeaf`](../../src/bartcore/model.hpp)'s calibration in R - a duplicated formula inside the check
  whose value is self-consistency.
- **C, a third entry, `sampleVarianceForestFromPrior`.** The leaf half is that class's own prior draw; the tree half is not
  one call. It needs `sampleTreesFromPrior`'s per-tree rejection draw under the empty-leaf veto and its attempt cap, then a
  `refreshVarianceForest` rebuilding `factorByTree` and `combinedVariance`, since
  [`Chain::setState`](../../src/bartcore/chain.hpp)'s variance validation demands every leaf a positive scale and every
  bottom occupied. No stream moves and both contracts stay true.

RECOMMEND C, priced at that tree half rather than as one entry. It does NOT lift
[`samplePriorPredictive`](../../R/dbarts.R)'s heteroscedastic `type = "ppd"` refusal: that function harvests through
`predict`, a mean-fit accessor, and the surface is readable only through a run's `variance` and `varianceTest` channels, so
lifting it needs a variance-surface accessor priced separately (landed: the accessor was folded into this slice by the
ruling below and the refusal went with it). The arm is the gaussian arm plus a variance-forest prior
draw per replication, ranking `s(x*)` off `varianceTest`, a VARIANCE on the original scale the functional square-roots
(`s.test` is the bart-object name, not the run channel's). A heteroscedastic aft arm follows once both land, reaching the
three-block cycle that [6. Gates](aft-variance-forest.md#6-gates) records as the honest gap. Unvalidated after all three:
the R-level survival packaging, left and interval censoring (unshipped, [Out of scope (v1)](survival.md#out-of-scope-v1)),
the flat header's log-likelihood channel, and - SBC being statistical at finite R - anything under its band.

## 7. Gates and poison

Nothing advances a generator when the target status is all events: the `GetRNGstate`/`PutRNGstate` bracket is an empty round
trip, `GaussianResponse::setResponse` ignores its `rng`, no sigma is drawn, and `refreshLatents` is guarded on a non-empty
censored set - so bitwise parity is reachable.

**(a) Creation parity, bitwise, at `updateScale = FALSE`.** A sampler created with an ALL-EVENT status must be bit-identical
to one created with a censored status and then set to all events at the same y. Pin `updateScale = FALSE`: at TRUE the scale
round trip multiplies and divides by `range_ * range_`, two roundings creation never performs, so the last ulp is
undefendable. Poison: let the setter rebuild the bounds but keep the old `censoredIndices_`; the redraw that should not
happen consumes draws, and the two runs part on sweep 1.

**(b) Ordering and bound exactness, no RNG.** In tests/cpp beside [`testAFTStateRoundTrip`](../../tests/cpp/test_model.cpp):
build A at `(y, S2)` and B at `(y, S1)`, set B to `S2`, and assert B's latents equal `y` at every event row, B's
[`AFTResponse::computeLogLikelihood`](../../src/bartcore/model.hpp) equals A's at every censored row - that entry reads the
bound, not the latent - and `sigmaDegreesOfFreedomForTesting()` did not move. Poison: run the status rebuild AFTER
`setResponse`, so a row censored under both takes its bound off its own latent; and count events rather than rows in the
degrees of freedom.

**(c) y and status in one call.** The only gate where both move: from `(y1, S1)` call `setResponse(y2, status = S2)`, then
assert every S2-censored row's bound is `y2` (by (b)'s log-likelihood comparison against a sampler created at `(y2, S2)`
over the same transform) and every S2-event row's latent is `y2`. Poison: build the bounds before the memcpy and they come
back at `y1`.

**(d) The redraw's law, and the R surface.** Extend [`testAFTCensoredMoments`](../../tests/cpp/test_model.cpp): after a flip
that newly censors a row, its redraws match the lower-truncated normal at ITS bound; poison, the pre-flip bound. In
[test-aft.R](../../inst/tinytest/test-aft.R): (a) through the handle; the refusals (off aft, a wrong length, a non-real
vector on the handle path, a value neither 0 nor 1, `NA`); and, R5 only, the handshake - set a status, store state,
invalidate the pointer, check the re-created sampler continues from the observed times, not stale latents. Poison: drop the
`survival.digest` write, and a row that was an event at store time and is censored after comes back sitting exactly at its
bound instead of being redrawn. That is the digest's own case, the reverse flip: the index-restricted `restoreLatents` is
what protects a row censored at store time and an event after, regardless of the digest.

## 8. Slices

1. **The setter.** [`AFTResponse`](../../src/bartcore/model.hpp) plus the new non-pure virtual; `Chain` and `Sampler`
   forwarding; the widened `bartcore_setResponse` and its `setData` message restatement, with the prototype in
   src/R_interface_bartcore.hpp and the arity in src/R_interface.cpp's `DEF_FUNC` table; the index-restricted
   `restoreLatents`, the `survival.digest` attribute and its reconciliation; the R5 method and its Rd; the handle call site
   [`bartcoreSetResponse`](../../inst/common/bartcoreHandle.R); gates (a) through (d). src/C_interface.cpp is untouched
   under (ii), edited under (i). A new `ResponseModel` virtual is a full-recompile hazard: `--preclean`. Roughly 260 lines
   of code and 190 of tests.
2. **The aft SBC arm.** Harness only: the `sbcFamilySpec` aft branch and `sbcFamilyConfig` entry, an entry in
   [`sbcMatrixConfigs`](../../benchmarks/R/sbc.R), the raised functional count with the M = 30 prose it stales, the measured
   burn and the workflow matrix row.
3. **The variance-forest prior draw and the heteroscedastic gaussian arm**, at the tree-half price.
4. **The heteroscedastic aft arm.** Harness only, once 2 and 3 land.

Matrix cells, in [feature-matrix.md](feature-matrix.md): the `setData` bullet's aft parenthetical and footnote f13 both
quote the message slice 1 restates; the Gaps row naming the setter as the SBC-coverage enabler closes with slices 1 and 2;
and the row calling heteroscedastic SBC coverage liftable via `setState` is OVERTURNED by slice 3, which the edit says
rather than restating.

**Landed.** Slice 1: fcd60feb (code, tests), e20c6462 (Rd, NEWS), f9bc9260 (null-status no-op); tip f9bc9260. Three
differences from the recommendation above: `AFTResponse` gained three facade virtuals rather than one -
`setSurvivalStatus`, `survivalDigest`, `reapplySurvivalStatus`, the `weightsDigest`/`reapplyWeights` shape; the `setData`
refusal keeps the fragment (["fix the censoring structure at creation"](../../src/R_interface_bartcore.cpp)) and restates
only its tail; and gate (d)'s [`testAFTCensoredMoments`](../../tests/cpp/test_model.cpp) gained an alone-path assertion,
since the joint call's own bound refresh hides the rebuild's bounds.

Slice 2's harness half: the reuse conversion, the [`sbcFamilySpec`](../../benchmarks/R/sbc.R) branch and
[`sbcFamilyConfig`](../../benchmarks/R/sbc.R) entry, the nine functionals and the NA burn. Two differences from section 6 as
proposed: the sampler is built through `dbarts()` with a `(time, status)` response, `dbartsSpec` being a specification
builder rather than a route to a sampler; and an empty censored set costs more than a dropped rank -
[`runSbcFamily`](../../benchmarks/R/sbc.R), [`rankUniformity`](../../benchmarks/R/sbc.R) and
[`sbcReport`](../../benchmarks/R/sbc.R) all had to learn to carry an NA rank, which is what "reported separately" is made
of.

Slice 2's admission. The ladder read 40000 sweeps over 24 prior-drawn datasets. Every functional clears ACF 0.1 by lag 39 -
sigma 39, `avg.f` 34, the five `f.star` cells 24 to 39, `S.star1` 35 and `logT.cens` 26, worst-case over the 24, with no
dataset leaving one undecorrelated - and the block-mean z is flat from the first block at that resolution; re-read in
400-sweep blocks only block 1 carries an offset (sigma +1.6, `S.star1` -4.9 mean signed z) and blocks 2 through 10 are flat.
So thin 40 covers the worst lag and [`sbcBurnSweeps`](../../benchmarks/R/sbc.R) takes 4000 sweeps, a 10x margin on a
transient that lives under 400 to 800.

At `R = 200`, `L = 150`, thin 40 all nine functionals PASS, and at the per-functional 5% band (0.0924), which is stricter
than the Bonferroni'd one the arm is admitted under: sigma 0.0620, `avg.f` 0.0635, the `f.star` cells 0.0377 to 0.0905,
`S.star1` 0.0736, `logT.cens` 0.0571. No replication drew an empty censored set, so `logT.cens` carried all 200 ranks. The
arm is in [`sbcMatrixConfigs`](../../benchmarks/R/sbc.R), [`sbcMatrixFunctionals`](../../benchmarks/R/sbc.R) is 39 and the
gaussian arm replays its ranks against the wider band unchanged. The workflow row (["config: aft"](../../.github/workflows/sbc.yaml))
runs thin 40 at a 15-minute timeout: section 6 priced 5-10 minutes at thin 30 and the run measures 79 s at the thin 40 the
ladder requires, so the timeout is the gaussian row's rather than the ~3x rule's, that job's floor being the dependency
install and the package build.

**Slice 3's prior-draw entry landed.** `Chain::sampleVarianceForestFromPrior`, recommendation C built at the price
section 6 states: the structure is the CGM prior conditioned on carrying no empty leaf, drawn by whole-tree rejection
against the user weights under `sampleTreesFromPrior`'s attempt cap and its one-scan settlement of the empty
conditioning event, and the factors are [`ConstantVarianceLeaf`](../../src/bartcore/model.hpp)'s own prior draw. The
rebuild is `refreshVarianceForest`'s order and arithmetic, so the entry leaves live state: the state a drawn chain
reports restores through [`Chain::setState`](../../src/bartcore/chain.hpp)'s variance validation. Both contracts stay
true - the two forest entries are untouched and every equivalence baseline replays bitwise. One deviation from section 6:
the recursion `sampleTreesFromPrior` uses was split so the variance forest, which owns a `CGMTreePrior` but no forest
object, reaches it directly; the mean path's call sequence is unmoved. Gated by
[`testVarianceForestPriorDraw`](../../tests/cpp/test_state.cpp), whose leaf arm scores the reciprocal factor - exactly
`chisq(nu) / nu` after scaling, so its two moments pin the calibrated scale and the calibrated degrees of freedom with a
closed-form standard error - and by the variance section of
[test-heteroscedastic.R](../../inst/tinytest/test-heteroscedastic.R). Poisons run: doubling the drawn factor fails both
moment arms; suppressing the structure draw fails both structure arms.

**The arm reads `s(x)` with a current-state accessor: the ruling, and what landed.** Section 6 specifies the
heteroscedastic gaussian arm as the gaussian arm plus this prior draw per replication, but a generator needs the drawn
`s(x)` at the TRAIN rows to simulate `y0` and at the test rows for its functionals, and section 6 itself records that the
surface is readable only through a run's `variance` and `varianceTest` channels - which report the state AFTER a sweep,
not the prior draw. The gaussian arm's own precedent did not settle it: that arm draws `sigma` in R because no engine
entry draws it, whereas here the draw IS in the engine and only the READ was missing. The maintainer ruled on 2026-09-13
for the first of the three candidates - a current-state accessor over [`Chain::varianceFits`](../../src/bartcore/chain.hpp)
and [`Chain::varianceTestFits`](../../src/bartcore/chain.hpp), folded into this slice rather than priced separately, with
[`samplePriorPredictive`](../../R/dbarts.R)'s heteroscedastic `type = "ppd"` refusal lifted in the same change. NOT taken:
routing the drawn trees in R off the reported state, which is option B's duplication moved from the calibration to the
surface; and deferring the arm until the accessor lands on its own.

**What landed.** [`dbartsSampler$getVariance`](../../man/dbartsSampler-class.Rd), `getVariance(test = FALSE)`: the
current variance surface `s^2(x)` on the ORIGINAL response scale, an n.observations x n.chains matrix at the default and
an n.test x n.chains matrix at `test = TRUE`. It reports exactly what a recorded sweep's `variance` and `varianceTest`
channels carry - the working product times `sigmaScale^2`, the scaling storeSample applies - at the same per-chain shape
the other current-state reads use ([`getFitsWithoutOffset`](../../R/dbarts.R)), and NULL exactly where those channels
report nothing: off a variance forest, and at a test read with no test rows. It addresses the trees IN FORCE rather than
the saved samples `predict` replays, which is what makes it answer after a prior draw with no `keepTrees`. The test arm
REBUILDS the test product before reporting it, that product being maintained only at a recorded sweep - so a
test-predictor swap or a prior draw cannot be reported stale. [`Chain::currentVarianceFits`](../../src/bartcore/chain.hpp)
owns both refusals, [`SamplerBase`](../../src/bartcore/facade.hpp) carries the virtual and
[`bartcore_getVariance`](../../src/R_interface_bartcore.cpp) is the one bridge entry; the shipped header gains nothing, a
state read being an R method under the pure-C rule
([pure-c-header](../plans/pure-c-header.md#pure-c-header)).

**The ppd lift.** [`samplePriorPredictive`](../../R/dbarts.R) at `type = "ppd"` on a heteroscedastic sampler draws the
variance forest from its prior beside the mean forest each sample and adds `s(x) eps` at the rows being predicted,
reading `s^2(x)` there through the accessor over the draw sampler's own test rows. `predict` cannot serve that read:
the private draw sampler forces `keepTrees` off and the drawn trees are never recorded. The homoscedastic path is
untouched and bit-identical, and no draw moves anywhere - all three equivalence baselines replay bitwise.

Gated by [`testCurrentVarianceRead`](../../tests/cpp/test_state.cpp), whose swap arm is the rebuild's own poison; the
facade coverage row ([`FacadeVirtual::currentVarianceFits`](../../tests/cpp/test_facade.cpp)); and the accessor section
of [test-heteroscedastic.R](../../inst/tinytest/test-heteroscedastic.R), which pins the channel agreement bitwise at one
and at three chains, both NULLs, the rebuild after a test-predictor swap, and the drawn factor's two moments against the
calibrated degrees of freedom and scale - the closed-form standard error
[`testVarianceForestPriorDraw`](../../tests/cpp/test_state.cpp) uses, read here through the accessor at a one-tree
near-zero-growth variance prior - plus the ppd's own per-column variance against the prior mean of `s^2(x)` the accessor
reports at those rows. ["ppd.variance"](../../inst/tinytest/test-prior-predictive.R) replaces the refusal's test.

**The SBC arms follow.** Slice 3's heteroscedastic gaussian arm and slice 4's heteroscedastic aft arm are harness-only
now, at section 6's shape with the generator reading the drawn `s(x)` here; nothing in either is blocked on an engine or
R-surface entry any more.

**Slices 3 and 4 landed: the two heteroscedastic SBC arms.** Harness only, at section 6's shape. `hetero` is the
gaussian arm plus a variance forest and `hetero-aft` is the aft arm with the surface where its shared sigma was;
neither ranks sigma, a variance forest pinning it, and both rank the surface in its place.
[`sbcConfigHetero`](../../benchmarks/R/sbc.R) and [`sbcConfigHeteroAft`](../../benchmarks/R/sbc.R) carry the designs,
[`sbcMakeHeteroSampler`](../../benchmarks/R/sbc.R) builds the one pinned sampler each reuses - a variance forest takes
`setResponse` only at `updateScale = FALSE`, which is the pin a reused arm wants anyway, and refuses `setSigma`, which
is why no branch calls it - and [`sbcHeteroPriorDraw`](../../benchmarks/R/sbc.R) is the prior draw the composition
needs: THREE entries, the two mean-forest ones being mean-forest ones by contract, reading the drawn f through
`predict` and the drawn `s^2` through [`dbartsSampler$getVariance`](../../man/dbartsSampler-class.Rd) at the training
rows, where it sets the simulated noise, and at the test rows, where it is a functional.

Functionals, [`sbcHeteroFunctionals`](../../benchmarks/R/sbc.R) carrying the two the arms share: `avg.f` and three
`f.star`, as the gaussian arm's; three `s.star`, the run's `varianceTest` channel square-rooted, that channel
reporting a VARIANCE; and `avg.log.s`, the mean of log s(x) over the training rows off the train-side `variance`
channel - a LOG because the leaf factors are multiplicative and their product is heavy-tailed, so the level of the
surface is ranked rather than the few rows carrying its tail. `hetero-aft` adds the aft arm's own two, both now read
at the row's own scale rather than a shared one: `S(t0 | x*)`, which is what the composition is for, and `logT.cens`.
Three test points rather than the gaussian arm's five: each carries two functionals here, and the matrix's band pays
for every one.

Wiring, [`sbcCheckVarianceChannel`](../../benchmarks/R/sbc.R) beside the other arms' checks. theta0's s comes from the
accessor while its posterior draws come from the run's channels, so accessor and channel must agree at one state: both
agree to 0, train and test. And the generator draws the mean forest before the variance forest, so what it reads of
the mean draw must not move when the variance draw follows: 0. `hetero-aft` runs
[`sbcCheckAftLatents`](../../benchmarks/R/sbc.R) as well, its fixture censoring 32 of 150 rows at the checked draw.

The ladders, 40000 sweeps over 24 prior-drawn datasets each, the aft arm's own resolution. Both arms read the same
transient: at 400-sweep blocks only `avg.log.s` carries a block-1 offset - mean signed z 1.56 on `hetero` and 1.88 on
`hetero-aft` over the 24 datasets, 0.38 and 0.69 at block 2, flat after - so both take 4000 sweeps of burn in
[`sbcBurnSweeps`](../../benchmarks/R/sbc.R), a 10x margin. Thinning is where they part. `hetero` clears ACF 0.1 by lag
39 worst-case over the 24 (`avg.f` 2, the three `f.star` 30 to 39, the three `s.star` 7 to 12, `avg.log.s` 22), with
no dataset leaving any functional above lag 40, so thin 40. `hetero-aft` reaches lag 97 (`f.star3`, on 2 of the 24;
`avg.f` 45 and `f.star1` 49 on 1 each, the three `s.star` 21 to 28, `avg.log.s` 38, `S.star1` 22, `logT.cens` 24), so
thin 100: censoring costs the MEAN surface roughly twice the lag at the same design, and costs the variance surface
nothing.

The verdicts, R = 200 and L = 150 at those thins. Every functional PASSES, and at the per-functional 5% band (0.0924),
which is stricter than the Bonferroni'd one (0.1347) the arms are admitted under. The ecdf statistic is what the verdict
reads; the rank histogram's chi-square sits beside it and is inside the matrix's alpha everywhere, `hetero`'s
`avg.log.s` closest at p 0.003.

| functional | `hetero` | `hetero-aft` |
|---|---|---|
| `avg.f` | 0.0479 | 0.0285 |
| `f.star1` | 0.0600 | 0.0559 |
| `f.star2` | 0.0641 | 0.0374 |
| `f.star3` | 0.0551 | 0.0344 |
| `s.star1` | 0.0494 | 0.0849 |
| `s.star2` | 0.0548 | 0.0541 |
| `s.star3` | 0.0711 | 0.0440 |
| `avg.log.s` | 0.0483 | 0.0427 |
| `S.star1` | - | 0.0595 |
| `logT.cens` | - | 0.0366 |

No replication drew an empty censored set, so `hetero-aft`'s `logT.cens` carried all 200 ranks. The runs measure 2.2
and 4.4 minutes single-threaded. Both arms are in [`sbcMatrixConfigs`](../../benchmarks/R/sbc.R),
[`sbcMatrixFunctionals`](../../benchmarks/R/sbc.R) is 57, and the workflow rows
(["config: hetero"](../../.github/workflows/sbc.yaml), ["config: hetero-aft"](../../.github/workflows/sbc.yaml)) take
the gaussian row's 15 minutes rather than the ~3x rule's, the aft row's reason: under a quarter hour the job's floor
is the dependency install and the package build.

The poisons, [`sbcHeteroPoisons`](../../benchmarks/R/sbc.R), opt-in through `SBC_POISON` as the latent BCF arms' are
and by-hand discrimination runs rather than recorded verdicts. "s-scale" simulates at TWICE the drawn s(x), the
surface's level wrong: on `hetero` it reddens every surface functional at once - `s.star` 0.8055, 0.7707 and 0.8919,
`avg.log.s` 0.9934 against the 0.1347 band - and leaves all four mean functionals inside it (0.0311 to 0.0870).
"s-shuffle" permutes the drawn scales across the training rows, so the level and the whole multiset are untouched and
only the row-to-row ASSIGNMENT is wrong - the per-observation channel the composition is made of: on `hetero-aft`
`avg.log.s` FLAGS at 0.2776 and `s.star3` at 0.1409, `s.star1` and `s.star2` press the band at 0.1004 and 0.1175 with
chi-square p 0.000 on all four, and every mean functional stays inside (0.0287 to 0.0884). The asymmetry is the
poison's own: theta0's `avg.log.s` is permutation-invariant, so what it catches is the posterior collapsing toward a
flat surface under scrambled data.

Both of those leave the PRIOR right and move the data away from theta0, so neither reaches the prior-draw entry
itself. "s-df" is the complement and the only one that does: the surface comes from a SECOND variance forest
calibrated at four times the arm's residual df while the fit keeps its own, and theta0 records exactly the surface the
data came from, so the single mismatch is the prior that surface was drawn from. On `hetero` it reddens the surface
and nothing else - `avg.log.s` 0.2118, `s.star2` 0.1717 and `s.star3` 0.1480 against the 0.1347 band, `s.star1`
pressing it at 0.1227, chi-square p 0.000 on the three and 0.004 on the fourth - with every mean functional inside
(0.0326 to 0.0647). It is what discriminates a miscalibrated
[`dbartsSampler$sampleVarianceForestFromPrior`](../../man/dbartsSampler-class.Rd), the one channel generator and fit
do not share.

A fourth mismatch is NOT offered, and its clean reading would be vacuous rather than reassuring: skipping only the
GENERATOR's variance prior draw, so theta0's surface is what the last replication's chain left. That REMOVES the
prior-draw entry from the generator rather than mismatching it. The surface it carries is marginally a prior draw only
if the run's posterior draws are already exact and the chain over replications has reached its stationary law - which
is what the arm exists to test - and a miscalibrated entry reads clean under it precisely because the generator never
calls the entry. It costs rank independence across replications besides, successive theta0 sharing a surface.

What closes. The honest gap [6. Gates](aft-variance-forest.md#6-gates) records - the joint calibration of (mean
forest, variance forest, censored latents), untested while aft was out of the matrix - is tested by `hetero-aft` and
closes. Section 6's list of what stays unvalidated after all three arms stands unchanged: the R-level survival
packaging, left and interval censoring, the flat header's log-likelihood channel, and, SBC being statistical at finite
R, anything under the band.
