# latent-subset-mask

agent: S0 sonnet (tests only); S1/S2/S3 opus (each adds an engine channel a
  family composes into its own precision vector, and S1 carries a bridge entry);
  S4 sonnet (Rd, NEWS, records, baseline re-record). Serialized: one
  implementer, each slice lands before the next starts.
rng: NEUTRAL on every shipped configuration, EVERY slice - nothing installs a
  mask today, and with none installed the family serves its pre-mask pointer BY
  IDENTITY, so `rollAndSetNodeAveragesFused`'s `forestWeights != nullptr`
  refusal (chain.hpp) is reached exactly as it is today. The trio therefore
  expects IDENTICAL at every slice and a deviation is a LEAK, never a re-record.
  DRAW-LAW-CHANGING on the opt-in masked path only, deliberately: a mask
  composes into `workingWeights()` and suppresses per-row latent draws, so a
  masked sampler consumes a different rng stream than an unmasked one. That is
  the feature. One shipped REPORTING change rides S1 (V3): a zero-composed-weight
  row's pointwise log-likelihood becomes NaN where gaussian wrote `-Inf`.
window: this arc's S1 must be preceded by the `empty-leaf-veto-fix` slice (see
  "Ordering against empty-leaf-veto-fix"), and the flat entry rides
  `docs/plans/dbarts-h-reshape.md` S1, whose preconditions are "gated not
  assumed" - V4, V5 and V8 are settled here, so reshape S1 may start. Pre-release
  by default. Does NOT reopen the empty-leaf veto's own six count-based sites,
  and absorbs no part of that fix.
budget: ~465 engine + ~95 bridge + ~50 R + ~160 docs code lines, ~135
  reshape-side lines carried by dbarts-h-reshape S1, and ~1310 test lines.
  Per slice: S0 ~140 test; S1 ~265 engine + ~95 bridge + ~50 R + ~560 test
  (tinytest ~330, tests/cpp ~230); S2 ~130 engine + ~290 test; S3 ~70 engine +
  ~180 test; S4 ~160 docs + ~140 test. Memory, unpriced because trivial:
  `numChains * n` extra doubles for the per-chain mask buffer and composite.
  Budgets are sized to the MANDATED ORACLES, not the engine delta.

(artifacts `.claude/latent-subset-mask-design/` - memo, critique, synthesis,
synthesis-verification. Gitignored; this plan is the tracked record and carries
the load-bearing facts rather than pointing at them.)

## Goal

A sampler accepts a per-observation 0/1 mask that says "row i is not in the data
set for this sampler this sweep", and honours it for the LATENT families zero
weights cannot reach. The row keeps its leaf occupancy and keeps receiving a
fit; it leaves every sufficient statistic, every family-level parameter update,
and its own latent draw. The channel is absolute and independent of
`setWeights`, it normalizes an all-ones mask to "no mask" in the ENGINE, and it
refuses anything that is not exactly 0 or 1. v1 covers gaussian, Student-t,
probit and ordinal.

## Context

The demand is the row subset that changes BETWEEN draws (principal
stratification, mediation, IV). `docs/design/r-c-division.md` "Demand headlines"
carries the corrected survey: princeBART builds six probit samplers and re-runs
`setData` on all six every outer iteration, five of them with a changed
latent-stratum subset, and it guards each block with `if (sum(<stratum>) > 0)`
so a stratum that empties is SKIPPED and that sampler silently keeps the
previous iteration's data. Its `dbarts_binary` helper widens dbarts's binary
detection to admit an all-0/all-1 stratum because a constant all-1 response is
not detected as binary today. Both are defects this channel removes by
construction. The same document's "Adoption slate" entry for the latent-family
subset mask is where that demand is priced, and it records the corrected
framing: Gaussian row subsetting ALREADY ships via zero weights
(`dbartsSampler-class.Rd`, the `weights` item: "Weights of zero exclude an
observation from the likelihood while keeping its fitted values"), so what this
arc opens is the latent families that refusal does not reach.

- The model claim, and its exact scope. The shipped refusals for probit
  (`R/spec.R`, `enforceBinaryWeightPolicy` in R_interface_bartcore.cpp) and
  ordinal are identification arguments about a weighted LIKELIHOOD: a
  truncated-normal latent scaled by `1/sqrt(w)` is not a coherent model, and
  `ProbitResponse`'s own class comment records that the pre-1.0 engine did
  exactly that and it was dropped. None of that argument reaches `w in {0, 1}`:
  a 0/1 weight is not a reweighted likelihood, it is the likelihood of the
  SUBSAMPLE, which both families have exactly at any row set. The design does
  not relax the refusal; it occupies the sub-case the refusal never covered.
  For LOGISTIC the claim is weaker - the shipped message rules on `w = 0` by
  name ("drop zero-count rows") - so logistic rests on the separate
  NON-REDUNDANCY argument in "Per family". For nbinom and aft the shipped
  refusals are not identification arguments at all (nbinom's is
  channel-routing, "exposure belongs in the offset"; aft's creation refusal
  gives no reason), so their slice rests on the derivation alone.
- The prior ruling, and its limit. `R/spec.R` extends an all-ones-is-absent
  courtesy to probit, ordinal AND nbinom, each setting `data@weights <- NULL`
  with the comment "weights identically 1 are the unweighted likelihood and are
  treated as absent". So the package has already ruled that a degenerate weight
  vector is a statement about MEMBERSHIP. The qualifier is load-bearing: the
  courtesy is CREATION-ONLY. Post-creation `refuseBinaryWeightChange`
  (R_interface_bartcore.cpp, reached from `bartcore_setWeights` and
  `bartcore_setData`) refuses every non-gaussian weight change REGARDLESS OF
  VALUE - measured: `$setWeights(rep(1, n))` on a probit sampler errors while
  `dbarts(..., weights = rep(1, n))` on the same family is accepted and
  normalized away. This design is what carries the ruling across the
  between-draws boundary; it does not inherit it there.
- The shipped courtesy is R-only. `R/spec.R`'s own comment says "the bridge
  keeps the same checks as a backstop". It does not: `enforceBinaryWeightPolicy`
  refuses ANY non-null weights for probit, all-ones included. That asymmetry is
  why this channel's normalizer goes in the ENGINE.
- `TResponse` is the composition precedent: it holds the borrowed
  `userWeights_` AND an owned `composite_`, writes `composite_[i] = w_i *
  lambda_i`, hands `composite_.data()` to the contained gaussian, and
  RECOMPOSES on `setWeights` and on `restoreLatents` (model.hpp). The mask takes
  exactly that shape.
- `bartcore_setForestWeights` (R_interface_bartcore.cpp) is the bridge template:
  capability probe FIRST, then length, then an O(n) validation loop, then
  forward. `Sampler::setForestWeights` (sampler.hpp) is the fan-out shape, and
  `SamplerShape::supportsForestWeights` (facade.hpp) is the
  derived-from-the-same-predicate shape for the capability bool.
- `Chain::weights_` is WRITE-ONLY after construction (five writes, no reader),
  so the response model really is the sole authority on served weights and the
  normalizer's identity restore is a one-line branch per family.
- Related plans: `docs/plans/zero-weight-exactness.md` (the shipped zero-weight
  semantics, the per-forest weight channel, and F1's re-association measurement),
  `docs/plans/sigma-df-zero-weights.md` (the observation-weight df precedent),
  `docs/plans/dbarts-h-reshape.md` S1 (which carries this arc's flat entry),
  `docs/plans/c-api-growth.md` (the reservation record and the return-polarity
  erratum), `docs/design/empty-leaf-veto.md` (the fix this arc waits on).

## Binding decisions inherited (do not reopen)

1. **The inactive latent draw is SKIPPED, not drawn-and-discarded** (except for
   Student-t; see rule 4). Forced, not chosen: the truncated-normal primitive
   (`ext_rng_simulateLowerTruncatedStandardNormal`) and Polya-Gamma
   (`ext_rng_simulatePolyaGamma`, src/external/random.c) are REJECTION samplers
   whose consumption depends on the bound and on psi, so drawing-and-discarding
   would desynchronize the stream and break T2(c) and T3a.
2. **Gaussian IS in scope.** T1 - `setActiveRows(a)` bit-identical to
   `setWeights(w * a)` - is the bitwise reference oracle the arc is built on.
   The comparator is `setWeights(w * a)`, NEVER "no weights".
3. **Multinomial masking is GLOBAL only.** `afterCombine` is data-free and the
   chain's `w` argument is DISCARDED for multinomial (combiner.hpp).
4. **The chain-level composition alternative is REFUTED** on the sigma-df leak
   (see "Why the family composes, not the chain"). Do not re-propose a
   chain-owned composite buffer.
5. **The mask does NOT ride the state block.** It follows `z` and the per-forest
   weight (`zero-weight-exactness.md` binding decision 3). Two sharp edges get
   documented instead; see "State carriage".
6. **This arc absorbs no part of `empty-leaf-veto-fix`** and touches none of its
   six count-based sites.

## What the channel is

`a_i in {0, 1}` per observation. `a_i = 0` means "row i is not in the data set
for this sampler this sweep". Not "row i has low precision", not "row i is
missing", not "row i is dropped from the design".

Mechanism, engine-side:

- `ResponseModel::setActiveRows(const double* a)` virtual, defaulting to
  `return false` (a refusal), overridden per family in scope - the shape
  `ResponseModel::setWeights`'s base no-op already uses.
- `Chain::setActiveRows(const double* a)` returns bool: it runs the single
  validating/normalizing O(n) pass, forwards to the response model, then
  invalidates the vector-leaf statistics cache.
- `SamplerBase::setActiveRows` / `Sampler::setActiveRows` fan out across chains
  and AND the results - the `Sampler::setForestWeights` shape (sampler.hpp).
- `SamplerShape::supportsActiveRows` is DERIVED FROM THE SAME PREDICATE the
  setter refuses on, so the probe and the setter cannot disagree; this is the
  `supportsForestWeights` precedent (facade.hpp, chain.hpp). A bridge that
  probes the CAPABILITY never switches on `shape.family` and never has to name
  a family - which is what lets Student-t work at all, since `shape.family`
  reports `gaussian` for it.
- Each family in scope owns `std::vector<double> activeRows_` (the raw mask,
  kept for recomposition) and a composite buffer, and returns the composite from
  `workingWeights()` while a mask is installed, or its pre-mask pointer (often
  `nullptr`) when it is not. **The null gate is POINTER IDENTITY**, so the fused
  path survives untouched whenever no mask is installed.
- Each `Chain` owns its own `ResponseModel`, so each chain owns its own mask
  buffer and composite. The `{0,1}` scan is value-only and therefore uniform
  across chains, which is what keeps the fan-out from landing half applied.

Surfaces:

- `dbartsSampler$setActiveRows(active)` beside `$setWeights` (R/dbarts.R),
  validating length n, no NA, all in `{0, 1}`; `NULL` clears.
- `dbarts:::bartcoreSetActiveRows` beside `bartcoreSetWeights` (R/bartcore.R);
  the bridge entry follows `bartcore_setForestWeights`'s shape exactly.
- One flat entry (see "The dbarts.h footprint").
- NO `bart()`/`bart2()` argument: a one-shot creation-time row set is already
  spelled `subset =`.

## Semantics of inactive: the seven uniform rules

1. The row contributes nothing to any leaf sufficient statistic, branch
   log-likelihood, birth-scan `sumWeights`, or leaf parameter draw. Mechanism:
   the family composes `a_i` into its own `workingWeights()`, so every consumer
   of that accessor inherits it (eighteen call sites in chain.hpp at writing).
2. The row STILL OCCUPIES its leaf. `numObservations()` counts it, the veto sees
   it (`logLikelihoodForBranch`, moves.hpp), the scan's `count` sees it,
   `collapseEmptyNodes` merges by count when weights are null (tree.hpp).
   Unchanged from today's zero-weight semantics; this is where
   `empty-leaf-veto-fix` interacts.
3. The row STILL RECEIVES A FIT. `totalFits` spans all n rows, so `run()$train`,
   `getForestFits` and `predict` report `f(x_i)` at an inactive row - the
   shipped zero-weight contract. This is what makes the channel worth more than
   compaction to the demand class.
4. The row's LATENT DRAW IS SKIPPED, **for probit, ordinal, logistic, nbinom,
   multinomial and aft**. No rng is consumed for it and `latents()` returns its
   last drawn value, which is STALE. Documented, not silent; the correct read at
   an inactive row is its fit, not its latent. Forced by binding decision 1.

   **EXCEPTION, Student-t.** Rule 4 does NOT apply to `TResponse`. Its per-row
   latent `lambda_i` is drawn unconditionally for every row; the mask
   ANNIHILATES THE VALUE through the composite, it does not suppress the draw.
   Three merits, in order. (i) The alignment is internal to the oracle: T1's two
   arms run on the SAME build, and the comparator arm draws lambda for every row
   (`TResponse::refreshLatents`'s draw carries no weight gate), so the masked arm
   draws too - this is not compatibility with any external or historical stream,
   since the zero-weight Student-t path has never been published. (ii)
   Reactivation freshness: an inactive row keeps a current fit (rule 3), so its
   annihilated lambda is drawn from a conditional evaluated at current state, and
   t therefore escapes the one-sweep staleness rule 5 names for the skipping
   families - there is no stale latent to reactivate onto. (iii) The cost is one
   gamma draw per inactive row per sweep, and it is SCALE-FREE work:
   `ext_rng_simulateGamma`'s consumption depends on `shape` and the drawn
   variates only - `scale` multiplies at the exits and never enters an
   acceptance test (src/external/random.c) - with `shape = 0.5 * (nu + 1)`
   row-invariant, so a zeroed rate moves lambda's value and nothing else.

   **Considered and declined:** make rule 4 uniform by changing the SHIPPED
   zero-weight t path to SKIP the lambda draw at `w = 0` in both arms. It is
   attractive - no exception clause, no wasted draws, and T1-for-t stays bitwise
   because both arms skip - and it is declined on two counts. It converts this
   arc from strictly RNG-NEUTRAL to draw-shifting on a SHIPPED configuration,
   which re-pins every hardcoded regression value and re-records every baseline
   covering zero-weight t, and forfeits the arc's cleanest verification property:
   bitwise inert until a mask is installed (T2(a)). And it reintroduces for t the
   UNBOUNDED-AGE latent staleness that annihilation removes - a row inactive for
   k sweeps would reactivate on a lambda k sweeps old.
5. **Reactivation is a one-sweep MODEL hazard, not just a read hazard.** The
   sweep runs trees first and `refreshLatents` after, so the first sweep after a
   row REACTIVATES conditions that sweep's tree moves and leaf draws on the row's
   STALE latent (probit/ordinal z) or on a lambda drawn while the row was out
   (t). It does not disturb the subsample posterior while the mask is held fixed,
   and nothing was found that makes it incorrect, but for the target consumer - a
   mask that moves every outer iteration - it is not transient. It goes in the Rd
   beside rule 4.
6. The row's POINTWISE LOG-LIKELIHOOD is NaN (V3).
7. The row drops from every family-level parameter update that is a sum over
   rows, enumerated per family below.

## Per family

**gaussian (v1).** `a_i` multiplies the case weight. No latent to skip. Exactly
`setWeights(w * a)` by construction, which is what makes it the bitwise
reference oracle (T1).

**INSTRUCTION, not outcome: `setActiveRows` MUST route the composite through the
counting path.** `numPositiveWeights_` is a cached scalar function of the
weights, set in `GaussianResponse`'s constructor and recomputed ONLY in
`setData` and `setWeights`, and `drawSigma` passes it into the posterior. A
gaussian `setActiveRows` that recomposes the buffer WITHOUT recounting
reintroduces exactly the sigma-df leak that kills the chain-level alternative -
inside the recommended design. Two acceptable implementations: call
`setWeights(composite_.data())` (which recounts), or recount explicitly in the
same pass. S1 pins it with `testSigmaPosteriorDf` (tests/cpp/test_model.cpp)
driven under a mask.

**Student-t (v1).** Recompose `composite_[i] = w_i * lambda_i * a_i` in
`refreshLatents` and in `recompose()`. The lambda draw CONTINUES for every row
(rule 4's exception).

The nu statistics must be CHANGED to read the COMPOSED weight. The gate inside
`TResponse::refreshLatents` reads the LOCAL `w`, which is the raw user weight,
not the composed one, so an inactive row with user weight 1 would stay in
`numInformative`, `sumLogLambda` and `sumLambda`, and therefore in
`ResidualDfPrior::drawIndex`. The gate becomes `w * a > 0`. The contained
gaussian inherits the df through the composite it is already handed, subject to
the same recount obligation.

**probit (v1).** Skip the truncated-normal draw in `refreshLatents`; the stale
`latents_[i]` stays finite, so `rebuildWorking` and the running residual never
see a NaN. `workingWeights()` returns the 0/1 vector while a mask is installed
and `nullptr` when it is not, preserving the fused path in the unmasked case.

**ordinal (v1).** As probit for `drawLatents`, PLUS two sums that are both
required: `computeScales` tallies only active rows - it changes the PROPOSAL -
and `cutpointLogAcceptance` sums only active rows - it changes the TARGET. The
MH accept/reject is per free cutpoint and reads only `gamma`, so it stays
stream-aligned with the compacted sampler; note that `updateCutpoints` consumes
1 or 2 variates per cutpoint, not always 2 (its out-of-order `continue`
precedes the accept draw).

**logistic (S2).** Skip the `PG(w_i, psi_i)` draw; leave `omega_[i]` and
`working_[i]` at their stale finite values; return `a_i * omega_i` from a
separate composite. **NEVER write the zero into `omega_`**: `working_[i] = w_i
(y_i - 0.5) / omega_i - offset_i` divides by it, and `0.0 * inf` in the node
kernels is the NaN hazard. The channel is NOT redundant with the count weights:
a zero count is refused at creation, no count change is accepted after it, and
`reps = 0` would still draw one PG variate unconditionally.

**nbinom (S2).** As logistic for omega, PLUS the dispersion block: `sumLog1mP`
sums active rows only and `NBDispersionPrior::computeKernel` is rebuilt over the
ACTIVE count histogram on every mask change. The `maxCount` scan must be
restricted to active rows too, or the grid length differs from the compacted arm
and T3a fails for a reason that is not a defect. That rebuild is
`O(n + maxCount * gridSize)` per install and is the one real per-install cost in
the design.

**aft (S2).** Compose into the contained gaussian's weights (so the sigma df is
inherited through the gaussian path) and skip the censored-row redraw for
inactive rows. Inherits T4(ii): the contained gaussian is constructed with
`nullptr` weights and rescales over all n rows.

**multinomial (S3, GLOBAL only).** Skip the K `PG(n_i, psi)` draws for an
inactive row in `drawForestGlue` and compose `a_i` into `forestWeights_` in
`formForestResponse` (combiner.hpp). Never write into `omega_`;
`formForestResponse` divides by it. `afterCombine` needs no change - it reads
leaf values and the prior precision only. PER-FOREST masking is REFUSED with a
model reason: the margin `C_if` is a log-sum-exp over the other K-1 forests, so
a row absent from category k's forest is still in every other category's
likelihood; "row i is out of category k only" restricts no likelihood.
Confirmed structurally by the chain's `w` argument being discarded for
multinomial, which is also why the mask cannot ride `workingWeights()` there.

**grouped (delegating).** `GroupedResponse` forwards `setActiveRows` to its base
exactly as it forwards `setWeights`; `drawGroupEffects` already weights
per-group sums by `w_i`, so an inactive row leaves its group's mean and
precision and an entirely-inactive group falls back to its prior through the
same formula. Needs a PIN, not an edit.

**heteroscedastic (decorating).** `formMeanWeights` (chain.hpp) divides the
composed weight by `s^2(x_i)`, so `0 / s^2 == 0`, and the variance forest's
suffstat carries the same composed weight. No edit; untested to date, so it
needs a PIN.

**BCF.** Gaussian response, so the mask composes into `GaussianResponse` for the
df and then multiplies into the per-forest weights at `composeForestWeights`.
A per-forest mask is refused as REDUNDANT: `setForestWeights` already expresses
it.

## Why the family composes, not the chain

The closest shipped precedent is `composeForestWeights` (chain.hpp), and
zero-weight-exactness S2 chose that shape for the per-forest weight. It is
refuted here on four counts, the first decisive:

1. **The sigma df leak.** `GaussianResponse` counts positive weights from its
   OWN pointer and `drawSigma` passes that count to the posterior. A chain-owned
   composite never reaches it, so a masked gaussian would draw sigma at the
   UNMASKED df: precisely the defect `docs/plans/sigma-df-zero-weights.md` fixed
   and `tests/cpp testSigmaPosteriorDf` pins at 4e5 draws. Family composition
   gets this right for free because the composite goes in through `setWeights`.
2. **Direct readers bypass a chain buffer.** `formMeanWeights`,
   `sweepVarianceForest` and SEVEN `collapseEmptyNodes` sites (including the one
   in `sampleTreesFromPrior`) call `response_->workingWeights()` directly. Each
   is a silent leak.
3. **It buys nothing where the work is.** The latent-draw skip, the ordinal
   cutpoint sums, nbinom's dispersion statistic and the grouped per-group sums
   are family-internal by nature.
4. `composeForestWeights` is reached only inside `if (combiner_)`, so on a
   single-forest chain the shipped machinery is not even on the path.

## Normalization and validation

ONE O(n) pass in `Chain::setActiveRows` does both jobs:

- `active == nullptr` -> clear, return true.
- any element not exactly `0.0` or `1.0` (NaN fails both tests) -> return false,
  install nothing.
- zero count == 0 -> **CLEAR** (drop any installed mask, restore the pre-mask
  pointer by IDENTITY), return true.
- otherwise -> copy into the family's owned mask buffer, recompose, invalidate.

The scan lives in the ENGINE, not in R and not only in the bridge, because only
the engine sits under all three surfaces. That is the lesson of r-c-division
defect 4, which is worse than a missing guard:
`dbarts_sampler_setWeights` (src/C_interface.cpp) has no family guard, no length
check, no finiteness check and no non-negativity check, and for probit/ordinal
it lands on `ResponseModel::setWeights`'s base no-op, whose comment ("the host
rejects earlier") is false for that one caller - so a flat caller's weights are
SILENTLY IGNORED. The O(n) cost does not cut against "fast over safe in C/C++":
`bartcore_setForestWeights` already runs an O(n) `R_FINITE && >= 0` loop in the
bridge before copying. R validates too (safe over fast in R, and it gets to
raise a good message); the bridge and the flat entry add only the capability
probe and the length check, then forward the bool.

Four consequences to document, none of them silent:

- `setActiveRows(rep(1, n))` reports SUCCESS and installs NOTHING; a mask that
  returns to all-ones CLEARS. `NULL` and `rep(1, n)` become the same object for
  masks exactly as they already are for weights.
- A fractional or NA element refuses the WHOLE call and installs nothing. No
  partial application.
- **All-zeros is ACCEPTED**, not refused. The zero count is n, not 0, so it takes
  the ordinary install path. Its behavior is coherent and is the documented limit
  of the shipped zero-weight semantics: gaussian sits at `numPositiveWeights_ ==
  0` (sigma drawn at the prior df), every leaf scores `logIntegratedLikelihood ==
  0.0` and draws from the prior, and probit/ordinal consume no latent rng at all.
  It is the natural answer for a host whose stratum has emptied - the very case
  princeBART guards around today - so refusing it would reintroduce the edge the
  design exists to remove.
- **The two channels carry OPPOSITE degenerate-value policies**:
  `setActiveRows(rep(1, n))` clears, while `setWeights(rep(1, n))` installs and
  measurably leaves the fused path. That asymmetry is deliberate - one channel is
  membership, the other is precision - and it gets its own Rd sentence rather
  than being left for a reader to discover.

## Cache invalidation

`Chain::setActiveRows` ends with the same two lines `Chain::setWeights` carries:

    if constexpr (L::hasVectorParams)
      forests_[0].leaf.invalidateStatistics();

on EVERY path that changes the installed mask, including the clearing path and
the normalizing path. This is required because the linear leaf's `U'WU` cache
re-validates on the ordered MEMBER LIST alone, so a mask toggle that moves no
membership is invisible to it, and `minCachedLeafSize = 32` means it bites on
the leaves that matter. `workingWeightsVaryPerSweep()` is false for gaussian,
probit, ordinal and aft, so the existing per-sweep invalidation does not cover
them; it is true for Student-t, where this addition is redundant but not wrong.
The leaf dispatch in `createSampler` (facade.hpp) has NO family gate, so
probit + linear leaf and ordinal + linear leaf are constructible.

`forests_[0]` alone is correct and is not a latent bug: both multi-forest
constructors static_assert `!L::hasVectorParams && !L::hasFunctionParams`, so a
vector-param chain has exactly one forest. The fix carries ZERO draw-law risk -
the cache's contract states a served value is bitwise the fresh scan's, so a
forced clear cannot move a draw - and is therefore unconditional rather than
conditional on the mask having changed.

**The GP leaf needs nothing.** Its zero-weight branches DO consult the same
cache, but every lookup re-compares the exact ordered member list
(`entry.members.size()` plus a `memcmp` against `tree.indices + node.begin`) and
the cached entry is a pure function of membership order and the FIXED
lengthscales - `gatherLeafCovariates` and `buildKernel` read no weights. A
weight or mask change cannot stale it because no weight ever entered it.
`anyZeroWeight` is recomputed inline from the live pointer at every node
evaluation and is never memoized, so the routing decision tracks the mask
unconditionally, and the zero-weight branch is the DOCUMENTED intent (its
comment: such rows "have infinite noise variance, so they contribute no
likelihood and the leaf takes the positive-subset paths below"). v1 SUPPORTS
masked GP leaves and pins them - but see the pin's SHAPE in S1: the two GP paths
consume different numbers of variates per leaf (`2 * numObs` clean against
`numObs + numPos` on the positive-subset path), so installing a mask on a
GP-leaf sampler moves the rng stream at every leaf and the pin must be
STRUCTURAL or statistical, never an equality.

## Composition with setWeights

Both setters are ABSOLUTE and INDEPENDENT. The family holds the borrowed user
weight pointer and the owned mask; `workingWeights()` serves `w * a` while a
mask is installed, `w` (by identity) when it is not. `setWeights(w2)` after
`setActiveRows(a)` KEEPS `a` and recomposes; `setActiveRows(a2)` after
`setWeights(w2)` keeps `w2`. Order-independent by construction. This is
`TResponse`'s shipped shape. The identity restore holds even over a pre-existing
user weight vector, because `w[i] * 1.0 == w[i]` exactly in IEEE.

`setResponse` / `setOffset`: no interaction; the mask survives them. `setData`:
the mask is length-n and `setData` may change n, so `setData` CLEARS the mask
and says so - consistent with `setData` cold-initializing the latents anyway.

## Grouped responses

Delegation, no refusal. The one thing to state in the Rd: a group all of whose
rows are inactive draws its effect from the PRIOR through `drawGroupEffects`'s
own formula - which is the coherent answer, and is NOT what a compacted sampler
with that group deleted would do.

## State carriage, and both sharp edges

The mask does NOT ride the state block (binding decision 5). No
`stateFormatVersion` bump, no `statesAgree` field, no `tests/cpp/common.cpp` or
`inst/common/stateContinuation.R` comparator.

Two edges, both documented in the Rd:

(i) A pipeline that REBUILDS a sampler from a stored state silently drops the
mask and computes different draws while `statesAgree` reports agreement. (The
shipped per-forest weight already has this edge.)

(ii) The LATENTS DO ride the state (the `"latents"` slot,
R_interface_bartcore.cpp), and under rule 4 an inactive row's latent is STALE
for every skipping family (probit/ordinal z, logistic/nbinom omega, aft's
censored draw). So the state block carries mask-history-dependent values, and
two samplers at the same posterior state with different mask histories report
DISAGREEMENT. Student-t is affected differently, not exempt: its lambda is drawn
every sweep for every row, but an inactive row's lambda is drawn from a
conditional its own data no longer informs, so the stored value is still
mask-history-dependent.

## Refusal messages

The engine returns false; the bridge raises. The honesty rule: every message
names WHICH of the three reasons applies - UNBUILT, IDENTIFICATION or
REDUNDANCY.

- `"active-row masking is not implemented for %s samplers"` naming the family,
  for logistic, nbinom, aft and multinomial in v1 - the reason is UNBUILT, and
  the message names the slice it is scheduled in.
- Permanent, per-forest multinomial: `"a multinomial mask applies to every
  category: the softmax margin reads all K forests, so a row cannot leave one
  category's likelihood alone"` - IDENTIFICATION.
- Permanent, per-forest BCF: `"a per-forest mask on a treatment (BCF) sampler is
  a per-forest weight; use setForestWeights"` - REDUNDANCY.
- Fractional values: `"active rows must be exactly 0 or 1: a fractional value is
  a weighted likelihood, which the latent families have no coherent form for"` -
  IDENTIFICATION, and it is the shipped argument in `R/spec.R` and
  `ProbitResponse`'s comment, unchanged.

Why the rule matters here: nbinom's shipped refusal is a channel-routing
argument ("exposure belongs in the offset") and aft's creation refusal gives no
reason at all, so neither is an identification argument this design can lean on.

## Ordering against empty-leaf-veto-fix

**Land the veto fix FIRST.** Both arcs re-record the zero-weight baseline;
sequenced fix-then-mask that is ONE re-record and the mask's baseline is recorded
under the final law. The reverse order needs two re-records and a documentation
amendment.

The two arcs meet exactly once: today a leaf all of whose members are inactive is
legal, scores exactly `0.0` (`ConstantGaussianLeaf::logIntegratedLikelihood`) and
draws its parameter from the prior at posterior precision 0, which a compacted
sampler could never produce. After the fix such a branch is vetoed.

**Struck, do not carry it forward:** the claim that the fix makes masked
occupancy MATCH compaction's state space. Five of the six count-based sites do
not move with the veto - the birth scan sentinels on `left.count`,
`collapseEmptyNodes` merges by count when weights are null, and the chi-k gates
are count-based. The fix changes which branches are VETOED, not which are
CREATED. The scheduling recommendation survives on the one-re-record reason
alone.

## The exactness target, per family

- **T1, BITWISE, gaussian and Student-t, SCOPED TO `a` NOT IDENTICALLY ONE.**
  For any mask `a` containing at least one zero, `setActiveRows(a)` on a sampler
  carrying weights `w` produces draws bit-identical to `setWeights(w * a)` with
  no mask, on every recorded channel. The comparator is `setWeights(w * a)`,
  NEVER "no weights" - an unweighted gaussian sampler is NOT bitwise equal to one
  carrying `rep(1, n)` (measured twice on two independent builds and two
  fixtures: `max|diff|` 1.179612e-14 and 5.329071e-15 on train, 3.330669e-16 on
  sigma; the mechanism is `rollAndSetNodeAveragesFused` refusing the instant
  `forestWeights != nullptr`, and gaussian-unweighted, probit and ordinal are
  exactly the three configurations still on the fused path).
  **The scope clause is not decoration.** At `a = rep(1, n)` the normalizer
  clears, so T2(b) governs that point and demands the no-weights sampler, while
  an unscoped T1 would demand `setWeights(rep(1, n))`; those two differ at
  ~5e-15, so an unscoped pair of pre-registered bitwise oracles is
  SELF-FALSIFYING and ships red on day one. T1 holds unscoped everywhere else:
  with `w` non-null the composite is elementwise identical on both arms and both
  have declined the fused path; with `w = NULL` and a zero present, the installed
  mask buffer IS the weight vector. T1 covers Student-t because both arms draw
  lambda for every row on the same build.
- **T2(a), BITWISE, every family.** With no mask installed, every family's draws
  are bit-identical to today. The null gate.
- **T2(b), BITWISE, every family.** An all-ones mask is bit-identical to no mask.
  Holds BY CONSTRUCTION after the normalizer, on the R5, bridge and flat surfaces
  alike. **Honest statement of what it covers: the NORMALIZER, not the composite
  path.**
- **T2(c), BITWISE, every family, per family.** Conditional independence:
  substituting arbitrary in-support values for the INACTIVE rows' responses (the
  0/1 label for probit, the category for ordinal, the count for
  logistic/nbinom/multinomial) leaves every ACTIVE row's recorded draw
  bit-identical. This is `test-zero-weights.R`'s existing oracle generalized, and
  it is the falsifier that PINS the skip semantics: because the latent primitives
  are rejection samplers, T2(c) FAILS if inactive rows' latents are drawn and
  discarded.
- **T3a, BITWISE, per family, in `tests/cpp`.** With identical inputs (eta,
  gamma, K, sigma held fixed), the masked family kernel over n rows and the same
  kernel over the compacted active rows produce bit-identical outputs at active
  rows and consume an identical rng stream. Covers probit `refreshLatents`,
  ordinal `computeScales` / `cutpointLogAcceptance` / `drawLatents`, logistic and
  nbinom omega, nbinom `sumLog1mP` and `computeKernel`, aft's censored redraw. It
  is bitwise because these loops are scalar, row-ordered and skip-shaped - they
  never multiply by a weight, so masking introduces no `0 * Inf`, no signed zero
  and no `log(0)` - and because none of them enters a `misc_` reduction. It also
  makes the ordinal comparison well-posed, since K is an input. **Four
  qualifications:** (i) T3a carries the model claim only TOGETHER WITH T2(a) and
  T2(c), because it cannot see the chain-level composition half - that `a`
  reaches every `workingWeights()` consumer - which is exactly where the leak
  risk lives; (ii) **sigma is EXCLUDED by name** -
  `ChiSquaredScalePrior::drawSigmaSqFromPosterior` goes through
  `misc_compute[Weighted]SumOfSquaredResiduals` (src/misc/moments.c), which is
  scalar but 5x UNROLLED and therefore re-associates with n, so no sigma-bearing
  statistic is ever bitwise between an n-row and an |A|-row arm; (iii) the aft
  arm requires the compacted comparator to PRESERVE ROW ORDER, since its censored
  redraw iterates gathered `censoredIndices_`; (iv) the ordinal kernels are
  PRIVATE and `tests/cpp` drives responses only through `createSampler`, so the
  ordinal arm needs engine-side `*ForTesting` hooks in the house idiom already
  used on that class (`cutpointLogAcceptanceForTesting`), priced in S1's engine
  budget. Gaussian and Student-t have NO T3a entry: they are carried by T1, which
  is what T1 is for.
- **T3b, STATISTICAL, ONE coarse arm per v1 family.** A masked sampler and a
  sampler freshly built on the retained rows over the SAME CUT GRID (a data
  handle view, R/bartcore.R) agree in the posterior mean of `f` at retained rows
  within Monte Carlo error. NOT bitwise, for three separately verified reasons:
  suffstat re-association (`zero-weight-exactness.md` F1, 1597/2000), count-based
  occupancy (rule 2), and a divergent tree-move rng stream. **Sized from a pilot
  written inside S1, not from this document.**
- **T4, NAMED DIFFERENCES, deliberately not closed.** Three, all in the Rd: (i)
  the cut grid is the FULL-data grid, not the active subset's - that is the
  FEATURE, it holds the split prior fixed while the row set moves; (ii) the
  response transform and the residual prior scale are the FULL-data calibration
  (`GaussianResponse::rescale` takes min/max over ALL n rows with weights
  ignored, and `sigmaSqPrior_.scale` derives from `range_`), so a masked
  gaussian - hence aft, hence any grouped or heteroscedastic decorator - keeps
  full-data `range_` and `sigmaSqPrior_.scale`; (iii) ordinal keeps the FULL-data
  K and its free cutpoints and log-gap prior even when a mask empties a boundary
  category, because `OrdinalResponse` fixes `numCategories_` at construction and
  nothing short of `setData` changes it. A consumer wanting the subset's own
  grid, calibration or K must still re-create.

## The dbarts.h footprint (carried by dbarts-h-reshape S1)

ONE entry, appended at the END of the X-list:

    X(int, dbarts_sampler_setActiveRows,
      (dbarts_sampler* sampler, const double* active), (sampler, active))

Contract (Doxygen):

- `active` is length `dbarts_sampler_numObservations`, each element exactly
  `0.0` or `1.0`; `NULL` clears (every row active); an all-ones vector is
  accepted and installs nothing.
- Returns **1 = accepted, 0 = refused** - the shipped convention
  (`dbarts_sampler_setPredictor` ends `accepted ? 1 : 0`; the polarity erratum in
  `docs/plans/c-api-growth.md`). No version constant moves (dbarts-h-reshape
  binding decision 8).
- **Ownership: the entry RETAINS NOTHING.** The values are consumed into the
  sampler's own buffer during the call and the caller's array is free
  immediately after it returns. This is NEITHER `setForestWeights`'s
  borrow-and-retain (which is why THAT entry obliges the caller to keep the array
  alive) NOR a copy into a holder; it is the predictor setters' "retain no
  pointer" clause. **No clause joins dbarts.h's keep-alive list.**
- Reachable for gaussian, Student-t, probit, logistic, aft, ordinal and nbinom -
  the families `dbarts_sampler_create` builds by name. multinomial and BCF have
  no flat creation path, so their masking is `dbarts:::`-only, as
  `bartcore_setCounts` and `bartcore_setForestWeights` already are.

Four obligations reshape S1 must carry:

1. **Validation layer.** The body runs the capability probe on
   `shape.supportsActiveRows` FIRST (never a family switch), then the length is
   implicit (there is no length argument - the entry reads `numObservations` from
   the shape), then it forwards the engine's bool. The exact-`{0,1}` scan is the
   ENGINE's, so the flat entry inherits it and the r-c-division defect-4 hole
   does not reopen. This is the only shape under which the `{0,1}` restriction is
   enforceable on the surface this arc ships.
2. **A body that ACCEPTS something.** The S1 body accepts GAUSSIAN (and
   Student-t, which is the same `shape.family`) and refuses every other family by
   name - the pattern reshape S1 PROPOSES at its item 5b: "The body accepts only
   what today's engine honours ... and refuses the rest naming the capability, so
   the family relaxes guard bodies later and moves no header." That item is
   itself CONDITIONAL on multiforest-extension-surface's fork 3, so it is a
   proposed pattern rather than shipped practice, and this obligation does not
   depend on 5b landing. Gaussian masking needs NO new engine work at that point:
   it is `setWeights(w * a)` composed at the entry.
3. **`consumer.c` + `test-capi.R` obligations.** Reshape S1's item 6 requires
   coverage of every entry plus the refusal matrix, and its gate list calls
   `test-capi.R` "the load-bearing gate". This entry owes: one positive arm
   (gaussian mask changes the fit, all-ones does not), one refusal arm per
   refused family reachable from `dbarts_sampler_create` (probit at S1 time), one
   fractional-value refusal, one `NULL`-clears arm. Because the body accepts
   gaussian, NO assertion has to invert when this arc's S1 lands - only refusal
   arms are relaxed, which is what item 5b was designed for.
4. **No dead entry.** A body that always returns 0 is WITHDRAWN: it would ship a
   symbol that lies by omission at the freeze and would plant an assertion in
   another arc's gate file that must later invert.

**Two entries considered and NOT proposed**, recorded so a later reader does not
reopen them: a reader (`dbarts_sampler_getActiveRows`) and a count
(`dbarts_sampler_numActiveObservations`). The channel does not ride the state
block, so the writer is the only source of the value and a reader can only echo
it. If V6 is ever reversed and the mask DOES ride the state, a reader becomes
necessary and must be added in the same re-bake.

## S0. Pins. No engine change.

`tests/cpp` cases asserting TODAY's behavior at the sites S1 moves:
`logIntegratedLikelihood` at `sumWeights == 0` returns exactly `0.0`; a leaf of
only zero-weight members is NOT vetoed; the scan's `count`/`sumWeights` split;
the GP leaf's `anyZeroWeight` routing. Plus tinytest pins that `$setWeights` on a
probit sampler still refuses and that `weights = rep(1, n)` is still normalized
away at creation.

Budget ~140 test. Gates: `tests/cpp` from clean. No `--preclean`, no trio (zero
engine delta). RNG class: NEUTRAL.

## S1. The channel + gaussian, Student-t, probit, ordinal

THE draw-law slice for masked configurations; bitwise-neutral for every unmasked
one. Carries the normalizer, the validation scan, the cache invalidation, the
composition order, the pointwise-log-likelihood decision, the R5 method, the
bridge entry, and T3a's kernel oracles.

Four instructions that are easy to read as outcomes and are not:

1. Route gaussian's composite through the `numPositiveWeights_` RECOUNT, and pin
   it with `testSigmaPosteriorDf` driven under a mask.
2. Change Student-t's nu gate to read the COMPOSED weight.
3. Add the ordinal `*ForTesting` hooks T3a needs (`computeScales`,
   `updateCutpoints`, `drawLatents`), beside the shipped
   `cutpointLogAcceptanceForTesting`.
4. Write the heteroscedastic, grouped and GP pins - the GP one STRUCTURAL or
   statistical, NEVER an equality, because a masked GP leaf's rng consumption
   changes per leaf.

Also in this slice: the all-zeros arm (finite output, forest at its prior, fits
still reported); the T3b pilot that sizes T3b's fixture; and a
`bench-sampler` arm, masked probit sweep vs unmasked probit sweep, on a quiet
machine - the number a princeBART-class consumer needs, and the honest cost of
leaving the fused path (the 1.41-1.54x fused-path figure is a cost of the
recommendation, not only an argument against alternatives).

Budget ~265 engine + ~95 bridge + ~50 R + ~560 test. Gates: the house battery
below, in full. RNG class: NEUTRAL on every shipped configuration - the trio
expects IDENTICAL and a deviation is a leak, ABORT - and DRAW-LAW-CHANGING on
the opt-in masked path. ASAN leg REQUIRED (new engine code, a new owned buffer
per family per chain, a new borrowed pointer). `test_shape.cpp` extended for the
new `supportsActiveRows` POD field.

## S2. logistic, nbinom, aft

Per "Per family". nbinom carries the dispersion-kernel rebuild and its
per-install cost note. The refusal messages for these three are relaxed from
UNBUILT to accepted.

Carried from the S1 independent review, inside this slice's test budget:

1. Strengthen the S1 heteroscedastic and grouped pins with the bitwise
   masked-vs-`setWeights(w * a)` oracle - both configurations are
   gaussian-reachable, so the cheap oracle applies. The shipped pins are
   finiteness-only and would not detect the mask failing to reach the
   variance-forest sufficient statistics or the per-group sums.
2. Add a sampler-level ordinal T2(c) arm beside the kernel-level T3a coverage
   (the kernel arm pins the skip semantics; the gap is sampler-level coverage,
   not correctness).

Budget ~130 engine + ~290 test. Gates: the house battery. RNG class: NEUTRAL on
every shipped configuration (trio IDENTICAL, deviation is a leak, ABORT);
draw-law-changing on the opt-in masked path. ASAN leg REQUIRED.

## S3. multinomial (global only)

Combiner-level. The per-forest refusal lands here with its model reason.

Budget ~70 engine + ~180 test. Gates: the house battery. RNG class: NEUTRAL on
every shipped configuration (trio IDENTICAL, deviation is a leak, ABORT);
draw-law-changing on the opt-in masked path. ASAN leg REQUIRED.

## S4. Surface, records, baselines

1. `man/dbartsSampler-class.Rd`, beside the `weights` item: the seven rules, all
   three T4 named differences, both state edges, the grouped all-inactive
   sentence, and the opposite-degenerate-policy sentence.
2. `inst/NEWS.Rd`, including the pointwise-log-likelihood REPORTING change.
3. A named recipe cross-referenced from `?bart` (feeding TODO `named-recipes`).
4. `docs/plans/c-api-growth.md`: the `dbarts_sampler_setActiveRows` entry with
   its retains-nothing ownership clause and its 1/0 polarity.
5. The equivalence-baseline re-record, plus two new `benchmarks/R/equivalence.R`
   scenarios (`maskprobit`, `maskordinal`, beside `zeroweights`) and one
   `bcf-equivalence` arm, ALL NULL-GATED so the EXISTING scenarios must stay
   BITWISE - any divergence there is a leak, ABORT.
6. The `## Landing` note in this file. There is no `docs/design/` document for
   this arc, so the plans README's Status-bump rule has no target; the record is
   this file plus the `docs/design/r-c-division.md` Adoption-slate entry.

Budget ~160 docs + ~140 test. Gates: full tinytest; `air format --check .`;
`lintr` on touched R files; `R CMD check` (Rd and NEWS move); the trio. RNG
class: NEUTRAL. The new scenarios are additions with no prior baseline; the
pre-existing ones must reproduce bitwise inside the same compare.

## Verification (the house battery, every engine slice)

- `R CMD INSTALL --preclean` into a PRIVATE library - MANDATORY: the
  `SamplerShape` POD gains a field and `facade.hpp` virtuals move, and stale
  objects bus-error. Delete the `benchmarks/kernels` binaries after any header
  edit (no dependency tracking there); `tests/cpp` tracks headers via `-MMD`.
- `cd tests/cpp && make clean && make && ./test_bartcore` - all pass. PLAIN AND
  ASAN legs (`ASAN_OPTIONS=detect_container_overflow=0`) for S1, S2 and S3, each
  of which adds new engine code and a new reachable buffer.
- `tests/cpp/test_shape.cpp` extended for the `supportsActiveRows` POD field
  (its `CHECK_SHAPE_FIELD` list). S1 only; S2-S4 add no field.
- Full `tinytest::test_package("dbarts")` from the preclean install, FAILURES ==
  0, with NO snapshot regenerated at any slice. A forced regeneration is a signal
  the slice reached further than intended - STOP and report.
- **The trio, EVERY slice, EXPECTING IDENTICAL.** This arc is bitwise-inert until
  a mask is installed and no baseline scenario installs one, so the expected
  verdict at S0, S1, S2 and S3 is exact match on all three, and any divergence is
  a LEAK, never a re-record - ABORT the slice. Use the CURRENT names per
  `benchmarks/baselines/MANIFEST`, which is authoritative; re-read it at spawn
  rather than trusting this line. At writing they are
  `equivalence-a825263.rds` (35 scenarios, "identical draws (same RNG stream)"),
  `bcf-equivalence-a825263.rds` (11 scenarios x their channels) and
  `multinomial-equivalence-1027be5.rds` (10 scenarios x their channels). No
  max-|z| line anywhere. S4 re-records only after adding its two null-gated
  scenarios, and states the PARTITION.
- **THE TRIO IS NECESSARY, NOT SUFFICIENT**: no trio scenario installs a mask, so
  it proves only that nothing leaked into the unmasked path. T1, T2(c) and T3a
  are the real oracles.
- `air format --check .` and `lintr::lint` on every touched R file (both enforced
  by `.github/workflows/lint.yaml`) - S1 and S4. `R CMD check` at S4.
- rchk on the next scheduled run: S1 touches `src/R_interface_bartcore.cpp`.
- Speed: S1 records ONE `bench-sampler` arm (masked vs unmasked probit) on a
  quiet machine, never concurrent with other load. It is a recorded NUMBER, not a
  pass/fail gate; no timed arm of `bench-sampler.R` installs a mask.

Stop conditions per `docs/plans/README.md`, plus:

1. Any trio scenario not reporting exact match -> STOP, the null gate leaks.
   Never re-record to clear it.
2. Any tinytest snapshot forced to regenerate -> STOP.
3. T1 failing on any arm where its scope clause says it should hold -> STOP; the
   composition or the recount is wrong. Do not widen the scope clause to make it
   green.
4. T2(c) failing -> STOP; the skip semantics are not what the design says, and
   the rejection-sampler argument means the stream has desynchronized.

## Falsifiers (pre-registered)

- **F1 (S1).** T2(c) must be shown RED against a deliberate draw-and-discard
  implementation of the inactive latent. A check that has never been red is not a
  check, and this is the one that pins binding decision 1.
- **F2 (S1).** T1 must be shown RED against a composite that skips the
  `numPositiveWeights_` recount - the sigma df is the whole reason family
  composition beat the chain-level alternative.
- **F3 (S1).** The vector-leaf (linear leaf) arm must be shown RED with the
  `invalidateStatistics()` call removed, on a probit or ordinal sampler with a
  leaf of at least `minCachedLeafSize` members whose membership does not move
  across the mask install.
- **F4 (S1).** T2(b) must be exercised on ALL THREE surfaces - the R5 method, the
  `dbarts:::` bridge entry, and the flat entry - because the whole point of
  putting the normalizer in the engine is that an R-only normalizer leaves the
  other two false.
- **F5 (S1).** A fractional element and an NA element each refuse the WHOLE call
  and leave the sampler's subsequent draws BITWISE unchanged. A refusal that
  perturbs is a bug.
- **F6 (S1).** The all-zeros mask RUNS - finite output, every forest at its
  prior, fits still reported - rather than refusing.
- **F7 (S1).** The GP pin must be structural or statistical and must be shown to
  DETECT a broken routing; it must NOT be written as an equality, because a
  masked GP leaf's per-leaf variate count changes by construction.
- **F8 (S2).** nbinom's `maxCount` scan restricted to active rows: T3a must be
  shown RED when the scan runs over all n rows, since the grid length then
  differs from the compacted arm.

## Edge cases the tests must name

An all-ones mask (clears, on all three surfaces). An all-zeros mask (accepted,
runs). A mask returning to all-ones after a real mask (clears; the pre-mask
pointer is restored BY IDENTITY, so the fused path comes back). A fractional
value, an NA, a wrong length (all refuse, install nothing). `setWeights` before
and after `setActiveRows`, in both orders, with the same final composite.
`setData` under an installed mask (CLEARS). `setResponse` and `setOffset` under
an installed mask (mask SURVIVES). A row reactivating after k sweeps (rule 5's
one-sweep hazard, and for t the absence of one). A group all of whose rows are
inactive. A heteroscedastic sampler under a mask. A GP leaf under a mask. A
linear leaf under a mask, with a cached leaf of at least 32 members. A masked
sampler whose state is stored and restored into a fresh sampler (the mask is
dropped; `statesAgree` reports agreement - edge (i)). Two samplers at the same
posterior state with different mask histories (`statesAgree` reports
DISAGREEMENT via the latents - edge (ii)). A per-forest mask attempt on a
multinomial and on a BCF sampler (refused, each with its own reason). An ordinal
mask that empties a boundary category (K, cutpoints and log-gap prior all
STAND - T4(iii)).

## NEWS bullets (inst/NEWS.Rd, at S4)

- A sampler accepts a per-observation 0/1 active-row mask (`$setActiveRows`),
  which removes rows from the likelihood and from the latent draws of a
  gaussian, Student-t, probit or ordinal sampler between draws while keeping
  their fitted values and their leaf occupancy. An all-ones mask is treated as no
  mask; fractional values are refused.
- The mask extends to logistic, negative-binomial and AFT samplers.
- The mask extends to multinomial samplers, where it applies to every category at
  once.
- The pointwise log-likelihood channel now reports NaN at any row whose composed
  weight is zero, in every family. Previously a zero-weight gaussian row reported
  `-Inf` and a zero-weight probit or ordinal row reported a finite value for a
  row that is not in the model.

## Decision record

Adopted at orchestrator discretion under VD's 2026-08-12
proceed-at-discretion grant. V1's Student-t disposition was individually
reviewed by VD the same day.

- **V1. Scope of v1: (a) gaussian + Student-t + probit + ordinal.** t is in on
  the MERITS, not on a correctness constraint - the override is ~15 lines, t is
  the family where zero-weight subsetting demonstrably already works, and it
  contributes a second bitwise T1 arm.
- **V2. The all-ones normalizer: (a) in the ENGINE**, inside
  `Chain::setActiveRows`, sharing the O(n) pass with the `{0,1}` scan. It is the
  only site covering the R5, bridge and flat surfaces alike, and it makes T2(b)
  true by construction; its cost is the T1 scope clause and the two-channel
  policy asymmetry, both documented above.
- **V3. Pointwise log-likelihood at inactive rows: (a) NaN** at every row whose
  composed weight is zero, in every family - following the channel's own shipped
  NaN flag, loud where rule 3 forbids silence. Costs a REPORTING change on a
  shipped gaussian configuration and a NEWS line.
- **V4. The flat entry at reshape S1: (a) append with a body that ACCEPTS
  gaussian** (and Student-t) and refuses every other family by name. Nothing
  inverts later; `test-capi.R` gets a positive arm.
- **V5. Fractional values: (a) refuse anything not exactly `0.0` or `1.0`,
  checked in the ENGINE's single O(n) pass.** The `{0,1}` restriction IS the
  identification argument, and the engine is the only site under all three
  surfaces.
- **V6. Does the mask ride the state block: (a) NO.** Follows `z` and the
  per-forest weight; no `stateFormatVersion` bump, no `statesAgree` field, no
  comparator changes. Both sharp edges are documented instead.
- **V7. Order against `empty-leaf-veto-fix`: (a) veto fix FIRST.** One baseline
  re-record, and v1 documents the final occupancy law.
- **V8. Name: (a) `setActiveRows` / `dbarts_sampler_setActiveRows`.** `setSubset`
  collides with `dbartsData(subset =)`, which PHYSICALLY drops rows and rebuilds
  the cut grid; `setMask` collides with the categorical rule bitmasks. The TODO
  entry and the design artifacts say "mask"; that vocabulary erratum is noted
  here rather than chased through them.
- **V9. Masked GP leaves: (a) SUPPORT and pin.** The zero-weight branch is the
  documented intent and the kernel cache cannot go stale, because it is keyed on
  the ordered member list it re-`memcmp`s at every lookup and holds no weights.
  The pin's SHAPE is constrained: structural or statistical, never an equality.

## Open items

- **T3b is PILOTED at S1** (gaussian; the pilot is not a shipped test). Sizing,
  from n = 500, p = 5 Friedman, 50 trees, 70% retained, burn 500, 4 replicates:
  - The two arms MUST take DIFFERENT rng seeds. At the same seed, 2 of 5
    replicates kept the masked and the compacted chain BITWISE LOCKED for 400
    sweeps (identical varcount, sigma to 2e-15) - the move stream only diverges
    once a re-association crosses an MH threshold, so a same-seed T3b is vacuous
    on some seeds and a real comparison on others.
  - The oracle is a RATIO against a within-arm control, not an absolute
    tolerance: rmse(posterior mean, masked vs compacted, at retained rows)
    divided by rmse(masked vs a SECOND masked run at another seed). Measured
    0.96-1.10 at S = 1000 and 0.99-1.21 at S = 4000, i.e. the cross-arm gap is
    the chain-to-chain gap. Threshold 1.5 leaves ~35% headroom. Correlation of
    the two posterior means, a cheap secondary, was 0.995-0.998.
  - A per-row z against a batch-means Monte Carlo SE does NOT work at any
    affordable S: it still reports rms|z| ~2.9 and max|z| ~10-15 at S = 4000
    while the ratio above says the arms agree. The SE is what is wrong (BART's
    posterior-mean autocorrelation), not the sampler. Do not pre-register it.
  - Cost: 0.15 s per arm at S = 1000, so ~0.5 s per replicate, ~2 s per family.
    T3b is cheap; S = 1000 and 4 replicates is the recommendation.
- **Nothing in this design was built or prototyped.** The cache-invalidation
  requirement is a source-and-contract derivation plus a reachability argument,
  not an end-to-end reproduction of a wrong draw; the heteroscedastic and grouped
  no-edit claims are arguments from `formMeanWeights`, `sweepVarianceForest` and
  `drawGroupEffects`, and each needs its S1 pin. If either pin fails, the
  composition claim in "Per family" is what fails.
- **The R-side all-ones courtesy is creation-only and stays that way.** This arc
  does not relax `refuseBinaryWeightChange`; a caller wanting a between-draws row
  set on a latent family uses the new channel, not `$setWeights`.

## Landing notes

S0 LANDED dc11a805, 2026-08-12. Two of the five planned pins were written:
the scan's `count`/`sumWeights` split (a bin of only zero-weight members
carries positive count and zero sumWeights, and the count-based occupancy
gate scores it; `tests/cpp/test_scan.cpp`) and the probit `$setWeights`
unconditional post-creation refusal, all-ones included
(`inst/tinytest/test-active-rows-pins.R`, new file, which S1 extends).
Three were already covered and not duplicated: `logIntegratedLikelihood`
at `sumWeights == 0` returning `0.0` (`test_model.cpp`
`testIntegratedLikelihood`), the GP leaf's `anyZeroWeight` routing
(`test_model.cpp` `testGPLeafZeroWeights`), and `weights = rep(1, n)`
normalized away at creation for probit/ordinal/nbinom (`test-family.R`,
`test-ordinal.R`, `test-nbinom.R`). The planned "leaf of only zero-weight
members is NOT vetoed" pin was DROPPED: the empty-leaf-veto-fix (21fc29c3)
landed after this plan was written and now vetoes such leaves; its own
tests pin the new law, and S1's spec sites should be read under it.
Gates: tests/cpp from clean, full tinytest FAILURES == 0; tests-only
slice, independent re-run delegated to CI on the push - six-green
confirmed 2026-08-12 (sanitizers, cpp-tests, exact-gates, R-CMD-check,
lint, pkgdown).

S1 LANDED 6db22aee, 2026-08-12. The channel (engine normalizer,
validation scan, cache invalidation, composition order), the R5
`$setActiveRows` + bridge entry, gaussian/Student-t/probit/ordinal, the
T3a hooks, the all-zeros arm, and `test_shape`'s `supportsActiveRows`
entry. V3 shipped on BOTH paths: engine `computeLogLikelihood` and the
R-side `pointwiseLogLikelihood` that `extract(type = "loglik")` uses,
so the NaN reporting change is surface-wide (S4's NEWS bullet must say
so). Budgets: engine/bridge/R on or under; tests +741 (~32 percent
over, under the 1.5x stop, oracle-driven). F1-F7 each shown
red-then-green; F7's original finiteness GP pin did NOT detect the
disabled `anyZeroWeight` routing and was strengthened to
structural/statistical. Recorded deviations, all reviewer-ruled: the
UNBUILT refusal ships the plan's literal message without a slice name
(no plan references in shipped code); F4's bridge arm runs end-to-end
through `bartcoreRun` on a fresh `bartcoreSampler()` handle (the
helper creates, never aliases); T2(b)/refusal arms use an UNWEIGHTED
sampler (on a weighted one an all-ones mask composes to
elementwise-identical values and the arm is vacuous); T2(c) is bitwise
at ACTIVE rows exactly as specified - an inactive row's reported fit
wobbles ~5e-15 under its unused response because `finalizeTotalFits`
reconstructs from the working response (pre-existing mechanism,
verified untouched; pinned at 1e-12). T3b pilot findings are under
Open items. Bench arm (VD quiet-machine grant 2026-08-12, armab
alternating rounds, loadavg clean): masked/unmasked probit sweep
ratio 1.092x at 80 percent active (spread 1.0865-1.0954; n = 10000,
p = 10, 75 trees, 1 chain/thread), 1.024x at 50 percent active; a
two-point extrapolation puts the fixed weighted-path cost near 1.14x
as the inactive fraction goes to 0 and break-even near 60 percent
inactive (extrapolation, not measurement). Gates: implementer battery
green; INDEPENDENT reviewer battery green from scratch (preclean
private lib, tests/cpp plain+ASAN, tinytest 4222/0, trio identical on
all three baselines, air+lintr; F3 independently re-falsified); x86
box leg green (tests/cpp plain+ASAN, tinytest 4210/0 incl. test-simd
and the new pins, equivalence.R statistical OK; the bcf/multinomial
harnesses are bitwise-only and cross-host inapplicable - TODO
equivalence-harness-statistical-mode); CI six-green on the push.
Carried to S2 (from the independent review): see the S2 section.
Carried to S4: normalize the plan-label comments in the arc's test
files.

S2 LANDED <commit>, <date>.

S3 LANDED <commit>, <date>.

S4 LANDED <commit>, <date>.

ARC COMPLETE <date>.
