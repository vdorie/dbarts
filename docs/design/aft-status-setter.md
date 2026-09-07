# An AFT censoring-status setter, and the SBC arms it enables

Status: LANDED - slice 1 (section 8), 2026-09-07 (fcd60feb, e20c6462, f9bc9260); slices 2-4 remain PROPOSED.

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
- **updateScale.** The status re-anchors nothing and every refusal stands, so a heteroscedastic aft still takes only
  `updateScale = FALSE` ([`refuseVarianceForestScaleUpdate`](../../src/R_interface_bartcore.cpp)), both R entry points
  defaulting it FALSE. Not an invariant: the transform is over `logT_`, which equals the observed times only just after a `setResponse` -
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
[`dbarts_sampler_create`, `DBARTS_C_API_HASH`](../../inst/include/dbarts/dbarts.h), so the hash does not move and no
`LinkingTo` consumer recompiles. None wants it - stan4bart on bartcore maps an `"aft"` token but fits gaussian and probit
only, calling no bridge entry, and bartCause on dbarts-1.0 has no survival code. The door is one X-macro entry shaped like
[`dbarts_sampler_setActiveRows`](../../inst/include/dbarts/dbarts.h), a capability-status `int` return refusing off aft;
adding it re-bakes the hash, pre-release a note rather than a migration. One consumer-visible message does change:
`bartcore_setData`'s aft refusal ("fix the censoring structure at creation") becomes false as stated. The refusal STAYS - a
whole-data replacement may change n, which the status is stated over - but the reason is restated to name the channel that
now serves the caller.

## 6. SBC

An aft arm already exists and REBUILDS its fit every replication ([`sbcMakeAftFit`](../../benchmarks/R/sbc.R) inside
[`runSbcAft`](../../benchmarks/R/sbc.R)), which is why it needs [`sbcAnchorScale`](../../benchmarks/R/sbc.R) to pin the leaf
prior and an offset to pin the shift against a transform every rebuild re-derives from `range(y0)`. It is in neither
[`sbcMatrixConfigs`](../../benchmarks/R/sbc.R) nor the workflow matrix
(["config: gaussian"](../../.github/workflows/sbc.yaml) is the nearest arm).

**The aft arm.** The setter converts it to the reused-sampler shape every other arm has: one sampler built once at a fixed
build response through `dbartsSpec(family = "aft")`; per replication a prior draw of `(f, sigma)`, then
`logT0 = f0 + sigma0 * eps`, `status0 = logT0 <= logC` against the censoring times
[`sbcAddCensoring`](../../benchmarks/R/sbc.R) pins, `y0 = pmin(logT0, logC)`, an overdispersed second prior draw,
`$setResponse(y0, status = status0, updateScale = FALSE)` and one `$run`. The transform never re-anchors, so both pins go. Functionals: `avg.f`, `sigma`, `f.star` at the five test points,
and two only this family has - `S(t0 | x*) = 1 - Phi((log t0 - f(x*)) / sigma)` at a fixed `t0`, the reported deliverable,
and `logT0[i]` at the LOWEST-INDEXED censored row, ranked against that row's posterior latents read with `getLatents` by the
per-sample `run(0, 1)` idiom the discrete functionals use. The censored set is random per replication and can be empty; such
a replication contributes no rank there and that functional's R is reported separately. It alone ranks the truncated-normal
imputation.

Budget: measure a burn ladder first, which needs an aft branch in [`sbcFamilySpec`](../../benchmarks/R/sbc.R) - draw, fit,
burnRun, sample - and not only an entry in [`sbcFamilyConfig`](../../benchmarks/R/sbc.R), since
[`sbcBurnLadder`](../../benchmarks/R/sbc.R) dispatches through the former. That branch IS the reuse conversion. A sweep is a
gaussian one plus a truncated-normal draw per censored row; at a gaussian-like cost and a t-like burn, R=200, L=150, thin=30
lands near 5-10 minutes, hence a 20-30 minute timeout at the workflow's ~3x rule. The arm's functionals raise
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
lifting it needs a variance-surface accessor priced separately. The arm is the gaussian arm plus a variance-forest prior
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
