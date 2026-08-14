# zero-weight-exactness

agent: S0-S2 opus (S1 moves a posterior; S2 adds an engine channel and a bridge
  entry). S3 sonnet (docs, records, surface tests). Serialized: one implementer,
  each slice lands before the next starts.
rng: S0 neutral (tests only). S1 POSTERIOR-CHANGING on the fixed-glue BCF
  configuration - the sampler stops targeting a 1e-9-perturbed model and targets
  the model the user wrote down - and shifting on every other BCF configuration
  (a sweep-0 initialization effect). No single-forest and no multinomial chain
  instantiates the touched combiner, so both of those baselines stay bitwise and
  a deviation there is a LEAK, not a re-record. S2/S3 neutral: the weight channel
  is null-gated and no shipped configuration installs one.
window: none. Everything is pre-release; the sister packages in
  /Users/vdorie/Repositories update in lockstep. Runs concurrently with
  bcf-public-surface, which owns the state-carriage question (see "Doors").
budget: S0 ~90 tests/cpp lines; S1 ~15 engine + ~200 test lines, plus a
  baseline commit; S2 ~110 engine + ~70 bridge + ~15 R + ~280 test lines;
  S3 docs + ~70 test lines.

## Goal

BCF's forest combiner stops faking a multiplier it does not have. A per-forest
multiplier that is zero (or indistinguishable from zero at R's own almost-equal
tolerance) produces an exact zero response and an exact zero weight instead of
being pushed up to +/- 1e-9 and divided by, so a row that carries no information
about a forest contributes exactly nothing to that forest's sufficient
statistics and the forest's reported surface at that row is the sum of its leaf
values with no cancellation. A caller-settable per-forest, per-observation
weight then exposes the same exclusion to callers who cannot express it through
the glue.

## The threshold: binding, VD 2026-08-10

VD's ruling: "consider the default tolerances from almost-equal type functions
in R and in Python and use those, to ensure that we conform to expectation."

The threshold is `sqrt(DBL_EPSILON)`, and it is a SNAP TO EXACTLY ZERO, not a
clamp. VERIFIED, not asserted (R 4.6.1 on this machine):

- `formals(all.equal.numeric)$tolerance` is `sqrt(.Machine$double.eps)`, value
  1.4901161193847656e-08. `testthat:::testthat_tolerance()` (edition 3) returns
  the same value; so does the default of `tinytest::expect_equal`, which is what
  this package's own suite compares with. `waldo::compare` defaults to
  `tolerance = NULL` (exact) and takes testthat's value when testthat calls it -
  do not credit waldo with the default.
- numpy's `numpy.isclose(a, b, rtol=1e-05, atol=1e-08)`: the absolute default is
  1e-08, within a factor of 1.49 (fetched from the numpy reference page for
  `numpy.isclose`; numpy is not installed here, so this is a documentation
  citation, not a measurement).
- `sqrt(DBL_EPSILON)` is `2^-26` EXACTLY, so the constant is representable with
  no rounding: write it `0x1p-26` with the decimal value in the comment.
- The package already uses this tolerance for an almost-equal test of its own:
  `validateSpec`'s proposal-probability check (R/spec.R) compares with
  `all.equal` at its default.

Why widening is licensed, and the honest limits of that claim. Three arguments,
in decreasing strength:

1. **The threshold is a condition-number cap, and `sqrt(eps)` is the unique
   choice that bounds the loss at half the mantissa.** The reparameterization
   forms `y' = r/m`, so `|y'| <= |r| / tol`. With the snap at `2^-26` the
   amplification is capped at `2^26`, and `finalizeTotalFits`' cancellation
   (`total = forestY - resid + mu`) and `rollTreeResidual`'s carry an absolute
   error of about `|r| * 2^-27`, i.e. `sqrt(eps)` relative - exactly half the
   available digits. Today's 1e-9 clamp permits amplification `~2^30` and an
   error `~1.2e-7` (the critique measured 4.3e-7 on tau fits of 0.35 after
   several sweeps). Widening the snap band therefore IMPROVES the worst-case
   reported-fit accuracy by about 15x, on top of making the `m == 0` case exact.
2. **Conformance.** A multiplier that R's `all.equal`, testthat and tinytest
   would all report as equal to zero is treated as zero. That is the ruling.
3. **The discarded information is bounded and small.** A row inside the band
   contributes at most `w * tol^2 = 2.2e-16 * w` of precision (below the ULP of
   any O(1) leaf weight total) and at most `w * tol * |r|` of weighted response.

HONEST LIMIT, recorded because the plan must not overclaim: the widening is a
monotone improvement only on `|m| < 1e-9` and at `m == 0`, where today's code
clamps. On `1e-9 <= |m| < 2^-26` today's arithmetic is unclamped, and the snap
replaces an exact leaf conditional with one that drops a contribution of order
`w * tol * |r|`. That band is precisely where the reported surface is already
unreliable (argument 1), the event has probability of order `2 * tol * density`
per `b` draw - about 3e-8 for a fixture-scale posterior - and the perturbation is
at the tolerance every comparison in the R ecosystem treats as zero. S1 gates
this rather than assuming it (see "Band-hit report").

Scope of the snap: it is a property of the REPARAMETERIZATION, not of the model.
`combinedFits` and `drawGlue` keep the exact `b0`; a snapped multiplier is not
written back into the glue. Say so in the comment.

## Tolerance-site classification (the sweep)

Every near-zero constant in `src/` classified. This document is the home of the
classification; do not restate it elsewhere.

| site | constant | class | verdict |
| --- | --- | --- | --- |
| `BCFForestCombiner::formForestResponse` (combiner.hpp) | 1e-9 | equality-with-division: the arc's subject | SNAP at `0x1p-26`, S1 |
| new per-forest weight channel (S2) | none | exact multiplicative factor, no division | NO tolerance; see the departure below |
| `monotoneTreeIsFeasible` (model.hpp) | `tol = 1e-9` | feasibility slack, ONE-SIDED, and LIVE | LEAVE; fix the stale doc comment. See "Open decision" |
| `monotoneIntegrate` (model.hpp) | `tol = 1e-12` | quadrature convergence criterion, not an equality test; halved down the recursion | OUT OF SCOPE |
| `sign * DBL_EPSILON` latent fallbacks (model.hpp, probit and ordinal) | DBL_EPSILON | sentinel value for a NaN draw, not a comparison | OUT OF SCOPE |
| GP kernel `nugget` (model.hpp) | 1e-6 | regularizer, not a comparison | OUT OF SCOPE |
| split-probability and proposal-probability sum-to-one checks (`R_interface_bartcore.cpp:1158`, `:1330`, `R/A_class.R:402`) | 1e-10 | IS an almost-equal predicate at unit scale, hand-picked | UNIFY at S3 (VD 2026-08-10) |

## Binding decisions inherited (do not reopen)

1. **`model.hpp:166` hardening**: `logIntegratedLikelihood`'s guard becomes
   `!(sumWeights > 0.0)`. Verified bitwise-neutral by the critique (the two
   predicates agree on `+0.0` and `-0.0` and differ only on NaN and negatives,
   both unreachable under the bridge's `>= 0` validation and `w m^2 >= 0`).
   Take it in S1.
2. **NO weight-aware empty-leaf veto.** Occupancy is count-based and load-bearing
   in six sites (`logLikelihoodForBranch` in moves.hpp, the scan's occupancy
   test, `afterCombine`'s L and M, the chi-k accumulation,
   `recoverVarianceLeafValuesBelow`, `mapCutPointsBelow`); making the veto
   weight-aware without making all six weight-aware splits the invariant
   `refreshVarianceForest` asserts. A zero-weight leaf scores exactly 0.0, a
   finite comparable number, so the veto (whose stated rationale is that
   `-HUGE_VAL` must beat any finite branch score) has nothing to protect
   against. A weightless leaf is a well-defined exact prior draw consuming the
   same single normal.
3. **State carriage DEFERS to bcf-public-surface**, which runs concurrently. The
   per-forest weight does not ride the state, exactly as `z` does not; no
   `stateFormatVersion` bump, no comparator field in `tests/cpp/common.cpp` or
   `inst/common/stateContinuation.R`. If that arc adds a BCF re-creation path,
   `z` and the per-forest weights become re-creation-visible together and are
   decided together, THEN.
4. **Do NOT change `BCFState`'s initial `b0 = 0.0`.** It is the honest neutral
   initialization and the exact-zero path is the correct evaluation at it.
5. **Do NOT edit `docs/plans/multiforest-veto-rate-falsifier.md`.** It contains
   no occurrence of "weight"; its M1 exactness is STRUCTURAL (column mask, not
   row weight) and is neither strengthened nor weakened here. The note that
   per-forest zero weights are NOT an opt-out from the empty-leaf veto lives in
   this document only.
6. **No `dbarts.h` symbol in this arc.** The flat API has no BCF creation entry
   (`dbarts_sampler_create` routes to `createHolder`, never `createBCFHolder`),
   so the setter would be unreachable. Reserve the name and signature in
   `docs/plans/c-api-growth.md` at S3.

## Open decision (VD)

**`monotoneTreeIsFeasible`'s `tol = 1e-9`: I recommend LEAVING IT, against the
commission's reading.** Evidence: its doc comment calls it "a component-test
predicate", and that comment is STALE. `MonotoneConstantGaussianLeaf::
drawFromPriorForTree` (model.hpp) calls it as the ACCEPTANCE PREDICATE of the
rejection sampler for the constrained prior leaf draw, reached from live engine
code at `Chain::sampleNodeParametersFromPrior` (chain.hpp), which the bridge
exposes as `bartcore_sampleNodeParametersFromPrior` and R reaches through
`samplePriorPredictive`. So changing that constant is a draw-law AND target
change on the monotone family, not a free ride. Two further reasons: the test is
one-sided (`mu < a - tol`), so widening it ADMITS leaf vectors that violate the
monotone cone by up to the tolerance - a loosened constraint, not an exactness
gain, and the opposite in spirit to a snap-to-zero; and its natural scale is the
leaf-value scale (`scale/k`), not unit, so an absolute machine-derived constant
is no less arbitrary there than 1e-9. ACTION TAKEN IN S1 REGARDLESS: correct the
stale comment to name the live caller. What would change the recommendation: a
decision that ecosystem conformance outweighs the loosening, in which case the
constant change needs its own slice, its own monotone gates (the monotone
tests/cpp cases plus tinytest), and an rng verdict of its own.

**The caller-supplied per-forest weight takes NO tolerance band** - a departure
from the commission's reading, stated so VD can overrule. A supplied weight is
multiplied in; nothing divides by it, so `s = 0.0` is exact by construction
(`w * m * m * 0.0 == 0.0` and the accumulated `w' y'` is exactly zero for finite
`y'`) and a small positive `s` is numerically benign. Snapping it would buy no
exactness at all (with `m != 0` the response stays `r/m`, so the reported-fit
benefit does NOT follow the weight channel), would introduce a cliff for a
caller annealing membership weights toward zero, and would put a per-row branch
on a path that needs none.

RESOLVED - VD ratified the no-band departure 2026-08-10, restating the
governing principle: constants are selected to conform to what a user expects,
and a deviation that reads as "normal" is fine. A supplied weight used exactly
as given is the normal expectation - no R modeling surface snaps a small
`weights=` entry to zero (`lm`, `glm` use it verbatim) - so here the BAND, not
its absence, would be the deviation.

## Corrected divergence shape (pre-registration; MEASURED, not predicted)

The critique built the fix in a scratch copy and ran the gates. Its numbers
replace the memo's prediction and become MANIFEST content:

```
default        MISMATCH in: mu, tau, glue, sigma, train      varcount IDENTICAL
restricted     MISMATCH in: mu, tau, glue, sigma, train      varcount IDENTICAL
glue_toggle    MISMATCH in: mu, tau, glue, sigma, train      varcount IDENTICAL
weighted       MISMATCH in: mu, tau, glue, sigma, train      varcount IDENTICAL
set_treatment  MISMATCH in: mu, tau, glue, sigma, train      varcount IDENTICAL
```

- **No tree structure moves anywhere.** `varcount` is bitwise identical in all
  five scenarios. The memo's "varcount divergence is the sharpest confirmation
  that tree structure moved" is FALSE and must not enter the record. The reason
  is that the MH ratio is `prior * transition * exp(delta)` with `delta ~ 1e-16`
  on a control-only branch, so the decision was already the prior's to within
  1e-16; the fix makes it exactly the prior's.
- The fixture's `varcount` channel is the MU forest only
  (`numVariableCountForests() == 1`, `variableCountForest(0) == reportedForest()
  == 0`, not overridden by the BCF combiner), so it could never have spoken
  about tau. Checked separately with `bartcoreForestVariableCounts(bc, 1L)`:
  mu AND tau variable counts are both bitwise identical old vs new. S1 commit 2
  closes that hole by adding the tau varcount channel to the fixture.
- **`glue_toggle`'s mechanism**: its `glue` channel differs in 1 of 3 entries,
  and that entry is `a`, not `b1`. Under `update.b = FALSE`, `b1` is never
  drawn, so the memo's "its `b1` draw reads tau fits that moved" is impossible.
  `a` moves because `drawGlue`'s a-block reads `y - b_z tau`.
- **Magnitudes** (max absolute divergence per channel over the recorded sweeps,
  measured): 1e-11 to 4.3e-6, against channel sds of order 0.25. The three
  exact-posterior oracles are UNMOVED to 4 decimal places. This is a hygiene and
  well-posedness fix, not a corrected bug; the MANIFEST must say so or a later
  reader will over-read "RE-RECORD REQUIRED, all five scenarios".
- **Mechanism clause for the MANIFEST**, corrected: `b0 = 0` holds at creation,
  AND after any `setState` or `installForests` whose donor carries `b0 = 0` -
  the glue DOES ride the state (`serializeGlue`; `statesAgree` compares
  `hasBCF, a, aVariance, b0, b1`).

## Context (seams, all read at 7759672)

- `BCFForestCombiner::formForestResponse` (combiner.hpp) is the subject; its doc
  comment ("|m_f| is floored to keep the division finite") states the guard's
  only purpose and must be rewritten to state the `m == 0` semantics.
  `forestMultiplier` returns `glue_.a` for forest 0 and `b1`/`b0` by `z` for
  forest 1.
- `forestResponse[i] = 0.0` is REQUIRED, not cosmetic: `finalizeTotalFits` and
  `rollTreeResidual` (chain.hpp) read `forestY[i]` arithmetically, and the node
  suffstat kernels accumulate `w * y[i]` where `0.0 * inf` is NaN.
- The ONLY two sites where a `ForestResponse` reaches a forest sweep are the
  sweep loop and `growForestFromRoot` (chain.hpp), confirmed by grep. In
  grow-from-root the weights are consumed by `Tree::setNodeAverages` and
  `growTreeFromRoot` immediately, so S2's chain-level multiply must be applied
  the moment `formForestResponse` returns, BEFORE the tree loop.
- `rollAndSetNodeAveragesFused` (chain.hpp) returns false the instant
  `forestWeights != nullptr`. A combiner always supplies weights, so BCF and
  multinomial never take the fused path; that is why both of their baselines
  stayed bitwise across the fused-suffstat re-record, and it is why S2 must
  refuse single-forest samplers (handing one a non-null pointer is a silent draw
  change plus a measured 1.41-1.54x slowdown on the most common configuration).
- `bartcore_setTreatment` (R_interface_bartcore.cpp) is the validation and
  ownership template: capability probe FIRST (never a forest count), then
  length, then copy into a holder-owned buffer, then hand the engine the
  borrowed `.data()`. `PROT_COUNT` is a fixed enum, so per-forest multiplicity
  cannot take a PROT slot; copy, retain nothing.
- The zero-safety census of the constant-leaf path is settled and its two
  contested hazards are REFUTED (upheld on re-derivation): `recoverVarianceLeaf
  ValuesBelow` performs no division, and `scan.hpp`'s
  `right = total - left` cannot go negative because `left` is literally an
  intermediate partial sum of the same left-fold over non-negative bin weights.
  The GP, linear and monotone leaf guards are unreachable from BCF
  (`BCFForestCombiner` static_asserts a constant leaf and `createBCFSampler`
  instantiates only `SamplerFacade<ConstantGaussianLeaf>`); keep them in mind as
  precedent, not as coverage.
- Sigma df: a per-forest weight must NOT reach `numPositiveWeights_`, and by
  construction it does not - per-forest weights live in the combiner's
  `ForestResponse` while the count lives on `GaussianResponse`, and `drawSigma`
  scores `combinedFits`. See `docs/plans/sigma-df-zero-weights.md` for the
  observation-weight precedent; do not restate its derivation. This is an
  assertion to falsify (F4), not a change to make.
- Chi-k: `chi-k-empty-leaf-count.md` gated the k accumulation on
  `numObservations() > 0` because a forced-zero empty leaf is not a draw from
  the k-scaled prior. A zero-weight leaf IS such a draw, so it belongs in
  `kNumLeaves` and `kSumSquaredParams`; the count-based gate already gives that.
  No change; the distinction earns one comment.
- Second consumer for S2, beyond the latent-treatment sensitivity class:
  zero-inflated log-linear BART keeps all n rows with augmented weights of zero
  (TODO, `multiforest-mutation-gaps` record: "the weight channel suffices").
- Design docs: `docs/design/bcf.md` (the glue and multiplier), `docs/design/
  forest-combiner.md`, `docs/design/model-space-survey.md` door 2. Related
  plans: `docs/plans/forest-combiner.md` (the mandated BCF gate set),
  `docs/plans/sigma-df-zero-weights.md`.
- Design artifacts (memo, blind refuting critique) sit at
  `.claude/zero-weight-exactness-design/`, which is GIT-IGNORED - this plan is
  the only tracked record, so it carries the load-bearing facts rather than
  pointing at them. Verdict STANDS WITH AMENDMENTS; all four blocking and all
  thirteen advisory amendments are adopted, with the departures recorded at the
  end of this file.

## S0. Pre-fix pin. No engine change.

`tests/cpp/test_sampler.cpp`, beside the existing BCF combiner cases:
`testBCFFloorContamination` drives `formForestResponse(1, ...)` with `b0 = 0`
and asserts TODAY's values - `forestWeights[control] == w * 1e-18` and
`forestResponse[control] == resid * 1e9`. A second numeric case pins the
control-only-leaf marginal: transcribing `ConstantGaussianLeaf::
logIntegratedLikelihood` at `m = 1e-9` with 100 unit-weight rows, `k = 2`,
`scale = 1/sqrt(25)`, `rv = 1` gives exactly 0.0 where the algebraic value is
1.43e-19, i.e. the two O(1) terms cancel to rounding. These assertions INVERT at
S1 (the `testCategoricalNeverSplit` precedent); their purpose is the pre-fix
evidence the MANIFEST cites.
Gates: `tests/cpp` from clean. No `--preclean`, no trio (zero engine delta).

## S1. The exact-zero snap. TWO commits. THE draw-law slice.

**Commit 1, engine.**

1. `formForestResponse`: a named constant `0x1p-26` (`sqrt(DBL_EPSILON)`,
   1.4901161193847656e-08) with the derivation in its comment, and a two-branch
   rule - `std::fabs(m) < tol` writes `forestResponse[i] = 0.0` and
   `forestWeights[i] = 0.0`; otherwise divide as today. Rewrite the doc comment
   to state the `m == 0` semantics ("this row carries no information about this
   forest"), the condition cap (`|y'| <= |r| / tol`, half-mantissa bound), and
   that the snap is local to the reparameterization - `combinedFits` and
   `drawGlue` keep the exact `b0`.
2. `logIntegratedLikelihood` guard to `!(sumWeights > 0.0)`.
3. Correct `monotoneTreeIsFeasible`'s stale "(a component-test predicate)"
   comment to name `drawFromPriorForTree` as its live caller. Comment only; the
   constant does not move.
4. Tests: invert S0's pins (both quantities exactly 0.0 at `b0 = 0`); pin the
   constant FROM BOTH SIDES with a band-edge pair - `m = 0x1p-27` snaps,
   `m = 0x1p-25` does NOT and yields exactly `resid / m` and `w * m * m`; assert
   `|forestResponse[i]| <= |resid| * 0x1p26` for every row on every path (the
   condition cap, and the standing `0 * inf` probe); and a control-only-leaf
   marginal that is now exactly 0.0 on both sides of a birth.
5. Test the 1.3(c) payoff by MEASUREMENT, not derivation (tinytest): fixed glue
   `(1, 0, 1)`, a binary moderator so every control row sharing a moderator
   value shares every tau leaf; assert their reported tau fits are BITWISE
   identical (spread exactly 0). On the old build that spread is 4.3e-7 over
   ~40 distinct values. State the scope honestly in the test comment: treated
   rows are NOT made exact, they keep the ordinary ~5e-16 cancellation.

**Band-hit report (mandatory, and the reason this slice is not neutral by
construction).** `a`, `b0` and `b1` are SAMPLED, so widening the threshold
widens the set of draws taking the exact path from `P ~ 1e-9` to `P ~ 3e-8` per
draw. Do not assume neutrality; measure it. (i) Committed screen: after the
compare, check every recorded `glue` channel of all five bcf-equivalence
scenarios for any `|a|, |b0|, |b1| < 0x1p-26` other than the `b0 = 0`
initialization. (ii) One-time verification, uncommitted: a temporary counter in
the snap branch counting hits with `m != 0.0`, run over the five scenarios;
report the count. Expected zero. A nonzero count is NOT an abort - it is a fact
that goes into the landing note and the MANIFEST entry, because it changes the
explanation of the divergence.

Gates, commit 1:
- `R CMD INSTALL --preclean` into a PRIVATE library (combiner.hpp and model.hpp
  are headers); delete the `benchmarks/kernels` binaries, which have no header
  dependency tracking.
- `tests/cpp` from clean, plain AND the ASAN leg
  (`ASAN_OPTIONS=detect_container_overflow=0`).
- Full tinytest. Expect the four BCF files (test-bcf.R, test-blocks.R,
  test-interactions.R, test-multi-forest-seam.R) to PASS UNCHANGED: a
  decimal-literal sweep over all four returns zero hits, so nothing pins a
  sampled value. A forced snapshot regeneration is a signal the change reached
  further than the combiner - stop and report.
- The trio, with these EXPECTED verdicts and ABORT clauses:
  `equivalence-c8f661a` 27/27 "identical draws (same RNG stream)" - any
  divergence means the change leaked out of the BCF combiner, ABORT the slice;
  `multinomial-equivalence-ec2a3d0` bitwise on 3 scenarios x 5 channels - any
  divergence, ABORT; `bcf-equivalence-99205ee` the exact five-scenario shape
  pre-registered above.
- **The mandated BCF oracle set** (`docs/plans/forest-combiner.md` step 4;
  `docs/plans/bcf-b-ridge.md`: these MUST STAY EXACT): `bcf-exact.R` quick AND
  full, `bcf-exact-restricted.R`, `bcf-exact-weak.R`. These are the only oracles
  in the arc that can distinguish "the draws changed" from "the draws are now
  wrong" - the bitwise trio cannot, by construction. They are maximally exposed:
  `bcf-exact.R`'s mode 1 pins the glue at `(1, 0, 1)`, i.e. `b0 == 0` on EVERY
  sweep, and mode 2a fixes `b` at `(0, 1)`. Measured on a patched build: all
  three PASS, gaps identical to 4 dp in all three modes (mode 1 E[mu] 0.0001,
  E[tau] 0.0003; 2a 0.0006 / 0.0004; 2b 0.0009 / 0.0001). Cite the mode-1 result
  in the MANIFEST; it is the positive evidence the re-record needs.
- No `bench-sampler` run: the change removes a division on a per-row path and
  adds a branch; it is not a hot-path shape change and no timed arm is BCF.

**Commit 2, baseline, separately.** Add the tau varcount channel
(`bartcoreForestVariableCounts(bc, 1L)`) to `benchmarks/R/bcf-equivalence.R` -
the re-record is the one moment adding a fixture channel is free, and the
fixture today has NO bitwise guard on the tau forest's tree structure. Then
`Rscript benchmarks/R/bcf-equivalence.R record benchmarks/baselines/
bcf-equivalence-<newhash>.rds`, plus a MANIFEST entry in the house format
carrying: the corrected mechanism clause; the per-scenario channel list with
`varcount` IDENTICAL; the magnitude band 1e-11 to 4.3e-6 and the "hygiene, not
bug" sentence; the `bcf-exact` mode-1 unmoved result; NEUTRALITY (at the same
commit `equivalence-c8f661a` is 27/27 identical and
`multinomial-equivalence-ec2a3d0` is bitwise on every channel); the band-hit
count; and the demotion of `99205ee` to `historical`. Verify the new baseline
compares clean against itself, and state the PARTITION for the added channel
(the five pre-existing channels are what the compare against `99205ee` covered).

## S2. The per-forest weight channel. Bitwise-neutral.

Semantics, stated once in the engine header and once in the R docstring: a
per-forest weight is a multiplicative PRECISION factor on forest f's own leaf
conditionals, composing with the observation weight -
`forestWeights[i] = w_i * m_f^2 * s_{f,i}`. It does NOT remove the row from
occupancy, from the empty-leaf veto, from the combination (the row still
receives `m_f f_f(x_i)`), or from the residual sigma df. `s = 0` means "this row
carries no information about forest f". **Two sharp edges must be documented, or
a consumer will be misled**: (i) with `s = 0` and `m != 0` only the WEIGHT is
zeroed - the response stays `r/m`, so S1's reported-fit exactness does NOT
follow the weight channel; (ii) the weight lives in the holder, not the state,
so a pipeline that REBUILDS the sampler silently drops it and computes different
draws while `statesAgree` reports agreement (a sharper edge than `z`, which is
required at creation).

1. Engine: `ForestCombiner::supportsForestWeights()` false by default, BCF
   override true (the `supportsResponseMutation` posture).
   `Chain::setForestWeights(f, w)` RETURNS BOOL: false when `combiner_ ==
   nullptr`, when `!combiner_->supportsForestWeights()`, or when
   `f >= forests_.size()`; otherwise stores a BORROWED pointer in a per-forest
   slot (nullptr clears) and returns true. `SamplerBase::setForestWeights`
   virtual and `Sampler` fan-out (the `setTreatment` pattern) propagate the
   bool. `SamplerShape::supportsForestWeights` is DERIVED FROM THE SAME
   PREDICATE so the shape and the setter can never disagree. The refusal is at
   the ENGINE, not only the bridge: `tests/cpp` constructs chains directly, and
   the reserved flat entry returns 1 on acceptance and 0 on refusal (ERRATUM
   2026-08-10 - this line originally said the inverse).
2. Engine application, chain-level, at BOTH `formForestResponse` sites and
   immediately after the call, before the tree loop: if and only if forest f has
   a weight installed, multiply into a chain-owned scratch buffer and repoint
   `forestWeights`. One implementation, one place to audit; it composes with any
   future combiner without a per-subclass override to forget. Cost: one O(n)
   pass and one length-n buffer per weighted forest per sweep, off the tree-move
   hot path. When nothing is installed the pointer passes through unchanged with
   no allocation, no copy and no arithmetic - that is the neutrality gate.
3. Bridge `bartcore_setForestWeights(ptr, forestIndex, weights)`: capability
   probe on `shape.supportsForestWeights` FIRST (never a forest count, so a
   multinomial cannot slip through), then `forest < numForests`, then
   `Rf_xlength == n`, then every element finite and `>= 0`. Ownership:
   `std::vector<std::vector<double>> ownedForestWeights` sized at creation - the
   NESTED shape is load-bearing, since a flat `numForests * n` buffer would
   dangle the engine's borrowed pointer on any resize while resizing the outer
   vector moves inner vectors without relocating their heap buffers. Copy the
   SEXP (the `ownedTreatment` precedent), do NOT `retain` it and do NOT take a
   PROT slot. Translate a false return to `Rf_error` (defense in depth: the
   probe should already have refused). Registry entry in `src/R_interface.cpp`.
4. R: `bartcoreSetForestWeights(bcSampler, forest, weights)` next to
   `bartcoreSetTreatment` in R/bartcore.R, `dbarts:::`-only. Validate R-side too
   (safe over fast in R).

Gates: `R CMD INSTALL --preclean` MANDATORY (the `SamplerShape` POD and
facade.hpp virtuals both move - stale objects bus-error); delete the
`benchmarks/kernels` binaries; `tests/cpp` from clean, plain AND ASAN (a new
heap buffer and a new borrowed pointer indexed by forest is exactly the
newly-reachable-engine-code case); `tests/cpp/test_shape.cpp` updated for the
new POD field; full tinytest with a new `inst/tinytest/test-forest-weights.R`
and NO snapshot regenerated; the trio BITWISE on ALL THREE baselines (bcf
against the NEW hash) - any divergence means the null gate leaks, ABORT;
`air format --check .`; `lintr::lint` on every touched R file (a fixed landing
step, three CI-red incidents); rchk on the next scheduled run.

## S3. Surface consistency, records, reservations.

1. tinytest: a `storeState`/`setState` round trip and an `installForests` leave
   an installed per-forest weight IN FORCE (it lives in the holder, not the
   state), with the comment saying that is deliberate and matches `z`. The
   CROSS-sampler case is DOCUMENTED, not asserted - say so explicitly, since the
   round trip only covers the same-holder case.
2. Docs: `docs/design/bcf.md` (the multiplier and glue section; and bump its
   `Status:` line per the plans README landing rule);
   `docs/design/forest-combiner.md`; `docs/design/model-space-survey.md` door 2
   (shape (2) SHIPPED, shapes (1) and (3) open and unscheduled, plus the
   non-interchangeability fact below); `TODO` under
   `multiforest-mutation-gaps`; `docs/plans/c-api-growth.md` (the reservation:
   `X(int, dbarts_sampler_setForestWeights, (dbarts_sampler*, size_t forest,
   const double* weights), (sampler, forest, weights))` appended at the END of
   the X-list (ERRATUM 2026-08-10: "nonzero on refusal" here was the INVERSE
   of the shipped convention - the entry returns 1 on acceptance, 0 on
   refusal, and the `DBARTS_C_API_MINOR` bump clause is void under VD's
   no-increment ruling; see c-api-growth.md); `weights == NULL`
   clears; per-forest weights must NOT be folded into a future BCF creation
   struct, because membership is resampled per sweep; and
   `dbarts_sampler_setWeights` must NOT be widened with a forest index, since
   the two channels differ on the sigma df); `inst/NEWS.Rd`; the Landing note in
   this file.
3. The recorded negative result, whose home is this file: **per-forest zero
   weights are NOT an opt-out from the empty-leaf veto**, because the veto counts
   observations. The veto-rate falsifier's M1 is about per-forest COLUMN masks
   and rests on a structural argument; the two documents use "mask" for
   different objects. Do not edit that document (binding decision 5).
4. **The one fact for the deferred shape fork, recorded not decided**: because
   occupancy is count-based, per-forest zero weight and physical compaction are
   NOT the same draw law even in the limit - under compaction the excluded rows
   do not exist, so their cut positions cannot be split on and their leaves
   cannot be created; under zero weight they can, at the prior. The three door-2
   shapes are therefore not interchangeable implementations of one semantic.
5. **Sum-to-one unification (VD 2026-08-10: unify).** The three hand-picked
   1e-10 almost-equal checks on user-supplied probability vectors -
   `R_interface_bartcore.cpp:1158` (rule proposal probabilities), `:1330`
   (split probabilities), `R/A_class.R:402` (validity) - move to the
   conventional tolerance. C-side: the same named `0x1p-26` constant as the
   combiner snap, sharing its derivation comment. R-side:
   `sqrt(.Machine$double.eps)` (all.equal's default), matching `R/spec.R`'s
   sibling check. Error messages unchanged. Tests: one tinytest per surface
   pinning the band from both sides - a vector mis-normalized by 1e-9 (refused
   today) is ACCEPTED, one off by 1e-7 is refused. This widens accepted INPUT,
   deliberately: anything newly admitted is input `all.equal` reports equal
   to 1.

Gates: full tinytest; `air format --check .`; `lintr::lint` on touched R files;
`R CMD INSTALL` refresh of the private lib for the bridge change (a `.cpp`, no
header, so `--preclean` is not required) and the trio BITWISE on all three
baselines - item 5 touches `src/`, and input validation must not move draws;
any divergence is a leak, ABORT.

## Falsifiers (pre-registered)

- **F0 (S1, STOP clause).** If any `bcf-equivalence` SCENARIO comes back
  bitwise-IDENTICAL at commit 1, the census is wrong: stop, do not re-record,
  re-derive who reaches the floor. Evaluated at the SCENARIO level as
  `bcf-equivalence.R compare` reports it (restated per CHANNEL it would misfire
  immediately - `varcount` is identical on all five) and on the FULL fixture
  only; the smallest divergence observed is 6.5e-12 (`glue_toggle` sigma) and
  `quick` halves the sweep count.
- **F1 (S2, replaces the memo's unachievable oracle).** With `s_i = 0` on a
  subset S, `node.sumWeights` and `node.sumWeightedResponse` must be BITWISE
  unchanged when the S rows' RESPONSES are replaced by arbitrary other finite
  values. Simulated at 0/2000 mismatches. Do NOT state it as
  "zeroed equals dropped": `misc_computeIndexedWeightedSufficientStatisticsFast`
  has a `length % 5` prologue and sums 5-wide groups, so removing rows
  re-associates the sum - 1597/2000 bitwise mismatches when simulated. If a
  complement-accumulation oracle is still wanted it must carry a tolerance and
  be labelled as such. Make the substituted values large and the same test
  doubles as the `0 * inf` probe.
- **F2 (S2).** `setForestWeights(f, rep(1, n))` produces a BITWISE identical run
  to not calling it at all - the null-gate and application-path oracle in one.
- **F3 (S2), the consumer test, BOTH halves.** Under the fixed glue
  `(1, 0, 1)`, a run with `setForestWeights(1, z)` is BITWISE identical to a run
  with no per-forest weight, on EVERY sweep INCLUDING sweep 0 - it does not wait
  for a glue draw, because after S1 the `b0 = 0` multiplier already writes
  exactly 0.0 in both channels and the chain-level multiply gives `0.0 * 0.0 ==
  +0.0` there and a bitwise-unchanged `w b1^2 * 1.0` for treated rows. The
  NEGATIVE half is mandatory: under the DEFAULT glue F3 must be expected to FAIL
  from sweep 1 on, because `b0` is then a.s. nonzero and the two routes exclude
  different information. Write both halves into the test or an implementer will
  generalize it and chase a phantom.
- **F4 (S2).** Zeroing a per-forest weight leaves `numPositiveWeights_` and the
  sigma posterior DF unchanged. Assert the df, not just the draw.
- **F5 (S2).** `setForestWeights(f, rep(0, n))` runs, produces finite output,
  and leaves forest f a pure prior sample. No NaN, no Inf, no refusal.
- **F6 (S2).** The refusal is a RETURNED FALSE at the engine level, asserted in
  `tests/cpp` on a single-forest sampler and on a multinomial sampler, and a
  refused call leaves the sampler's subsequent draws BITWISE unchanged - a
  refusal that perturbs is a bug.

## Edge cases the tests must name

All rows excluded from forest f (every leaf a prior draw; the forest still
contributes `m_f f_f(x_i)`; a weight channel can say "uninformed", never
"absent"; F5 pins that it runs rather than refusing, since a refusal on all-zero
would break a legitimate transient during a membership resample). Weight toggled
mid-chain (supported: the combiner caches nothing per-forest across sweeps; a
leaf non-empty under S can be weightless under S', and nothing needs to happen
because occupancy is unchanged - this is where the count-based veto earns its
keep, since no membership change can strand an empty leaf). A negative or NaN
per-forest weight (refused at the bridge; the engine adds no assert, matching
observation weights - fast over safe in C/C++). Forest index out of range
(bound-checked against `numForests`; the variance forest is a separate member
and unreachable from BCF anyway, `createBCFSampler` returning nullptr on
`numVarianceTrees > 0`). `n` change under a stale weight vector (impossible:
`refuseMultiForestMutation` refuses `setData` at `numForests >= 2`). Saved trees
/ `keepTrees` / predict (a weightless leaf's prior-drawn value is stored and
replayed like any other; no new branch). `samplePriorPredictive` /
`sampleTreesFromPrior` (unreachable from a BCF handle today, and correct as is:
`sampleTreesFromPrior` reads `response_->workingWeights()` and never calls
`formForestResponse` - a prior tree draw should be membership-blind; document
the decision, do not thread the weight through it).

## Doors held open (recorded, not scheduled)

- **The 1e-10 sum-to-one checks: RESOLVED into S3 item 5** (VD 2026-08-10:
  "Unify. The reasoning is why."). Held here originally because widening the
  tolerance changes which user INPUT is accepted; VD accepted that consequence -
  anything newly admitted is input `all.equal` reports equal to 1.
- **`monotoneTreeIsFeasible`'s tolerance**, per "Open decision".
- **The two `-Inf` sites at zero USER weight in the pointwise log-likelihood
  channel** (model.hpp, gaussian and one family decorator), unguarded and
  pre-existing. Off this path - BCF sets `logLikelihoodIsDefined() == false`.
  Log for the TODO; do not fix here.
- **The mutable per-forest row MASK and physical COMPACTION shapes**, their
  cut-grid question and the O(|S|) cost reduction. Unscheduled; this arc
  prejudges nothing about them beyond the non-interchangeability fact in S3.

## NEWS bullets (inst/NEWS.Rd, one per slice, same commit)

- S1: a BCF forest whose scale multiplier is zero now contributes exactly zero
  weight and exactly zero response for the affected rows instead of a
  near-zero floor, so the treatment surface reported at control rows under a
  fixed glue is exact rather than accurate to about seven digits.
- S2: a BCF sampler accepts a per-forest, per-observation weight, so a caller
  can exclude rows from one forest's leaf conditionals without excluding them
  from the model. `dbarts:::`-only; no public surface in this release.
- S3: probability vectors for the tree prior are validated with the standard
  `all.equal` tolerance instead of a hand-picked 1e-10.

## Departures from the memo and the critique (record)

1. **The floor is not lowered to `sqrt(DBL_MIN)` and the guard is not kept at
   1e-9**; it becomes a SNAP at `sqrt(DBL_EPSILON)` per VD 2026-08-10. This
   supersedes the memo's recommendation and the critique's alternative.
2. **The critique's A5 test is now wrong and is replaced.** It proposed pinning
   the lowered floor with a `m = 1e-12` case asserting the unfloored `resid / m`;
   under the snap, 1e-12 is INSIDE the band and snaps to zero. The band-edge pair
   (`0x1p-27` snaps, `0x1p-25` does not) pins the constant from both sides
   instead, which is strictly stronger. A5's denormal caveat lapses with the
   `sqrt(DBL_MIN)` floor it was written against.
3. **The commission's "monotone changes no draws" reading is refuted**, with
   evidence, in "Open decision": `monotoneTreeIsFeasible` is the acceptance
   predicate of the constrained prior's rejection sampler, reachable from R.
4. **The commission's "widening is a monotone improvement" reading is amended,
   not adopted verbatim.** It holds on `|m| < 1e-9` and at `m == 0`; on
   `[1e-9, 2^-26)` the snap replaces exact arithmetic with a bounded
   approximation. The stronger and correct justification is the condition-number
   cap (half-mantissa bound), which the memo and the critique both missed and
   which makes the WIDE threshold better than the narrow one for the reported
   surface by about 15x.
5. **No tolerance band on the caller-supplied per-forest weight**, against the
   commission's reading; reasoning in "Open decision".
6. **The claimed waldo default is not real** (`waldo::compare` defaults to
   `tolerance = NULL`, exact). Verified defaults are R's `all.equal`,
   `testthat:::testthat_tolerance()` and `tinytest::expect_equal`.
7. **A13 is adopted with a placement change**: the tau varcount channel is added
   in S1 commit 2 (the re-record), not deferred, and the MANIFEST states the
   partition for it.
8. **A8 is adopted in full** (do not edit the veto-rate falsifier doc). Its
   write-conflict rationale has lapsed - that measurement's freeze is committed
   at 7759672 - but its substantive rationale stands: nothing in that document
   is wrong, so there is nothing to correct.

## Landing notes

S0 LANDED 4979656; S1 LANDED c820227 (engine) + d53efe6 (baseline re-record,
`bcf-equivalence-c820227.rds`; recorded in-worktree under the pre-rebase hash
and renamed at the rebase onto the docs commits). Implementer gates green per
commit; the orchestrator independently re-ran the battery on the final tree:
install `--preclean` into a fresh private lib, `tests/cpp` from clean plain AND
ASAN (`detect_container_overflow=0`), tinytest 3753/0, `equivalence-c8f661a`
27/27 identical draws, `multinomial-equivalence-ec2a3d0` 3x5 identical,
BCF 5x7 identical against the new baseline, `bcf-exact.R` quick OK
(mode 2a/2b gaps at the plan's values), `air format --check .` clean, lintr
clean on touched files (one pre-existing brace-style hit in the harness,
outside the diff).

Band-hit report: count 0 over the five scenarios (temporary instrumented
build, reverted, restored engine byte-identical); the only in-band glue value
anywhere is `glue_toggle`'s fixed `b0 = 0`, so the re-record divergence is
entirely the `b0 = 0` initialization and none of it is the widening.
Divergence exactly as pre-registered: five scenarios mismatch in mu, tau,
glue, sigma, train; `varcount` identical in all five, tau varcounts measured
identical pre/post as well - no tree structure moved. Magnitudes 6.46e-12 to
4.34e-06. Payoff measured: control-cell tau fits bitwise identical post-fix
(spread exactly 0, vs 3.67e-7/3.95e-7 pre-fix over 67/78 distinct values);
treated rows keep the ordinary cancellation, as scoped.

Accepted deviations: S0's marginal pin asserts digit-loss plus a nonzero
birth delta rather than a literal 0.0 (exact cancellation is
fixture-dependent - measured -1.78e-15 against an algebraic 1.19e-17); the S0
case renamed `testBCFZeroMultiplierSnap` at S1; the MANIFEST partition names
SIX pre-existing channels - this plan's "five" was wrong, the 99205ee compare
did cover `varcount`; two extra pre-fix builds taken for evidence (payoff
spread, tau varcount cross-check). Pre-existing fact found, not this arc's:
`./test_bartcore sampler` run suite-filtered fails
`testBCFGrowForestFromRoot`'s pinned characteristic value on old and new
engines alike (the pin depends on the shared runif01 stream position); the
full run passes. `benchmarks/README.md` no longer names baseline hashes
(MANIFEST is authoritative).

S2 LANDED 153d1dd (rebased over the composition-probe records; recorded
in-worktree as 98c1817). All six CI workflows were green on the S1 push
before this landed. The slice survived a mid-gates session kill: a second
implementer audited the inherited uncommitted work against this plan and
found three gaps, all fixed before any gate ran - F4 asserted a
positive-count proxy rather than the sigma df (fixed with a
sigmaDegreesOfFreedomForTesting hook pinning nu_0 + #{w_i > 0} == 3 + n
across an all-zero install), the bridge lacked the Rf_isReal conjunct (an
integer vector reaching REAL() was UB), and the stale-weight-under-setData
edge case was untested. Gates double-run (implementer, then orchestrator
independently on the rebased tree): install --preclean, tests/cpp plain +
ASAN from clean, tinytest 3779/0, trio bitwise (27/27 + 3x5 + 5x7 vs
c820227 - the null gate holds), air 0, lintr 0 on both touched R files.
Falsifiers: F1 bitwise-invariant to 1e300 substituted responses on zeroed
rows; F2 all-ones identical to no call; F3 both halves (fixed glue
identical including sweep 0; default glue divergent after sweep 1); F4 the
df pinned directly; F5 all-zero finite with a live tau posterior; F6
engine-level refusals on single-forest and multinomial with
refused-call-perturbs-nothing. Budget: engine 127 raw/~95 code (~110
budget), bridge 47/~35 (~70), R 28/12 (~15), tests 465/337 (~280; the
overrun is mandated commentary, no test cut - the sized-to-oracle budget
lesson again).

S3 LANDED 8c04f6e (rebased over the multiforest plan commit; recorded
in-worktree as 5a52afd). Same-holder round trips (storeState/setState and
installForests) proven to leave an installed per-forest weight in force;
the cross-sampler drop documented, not asserted, as specified. Docs
updated: bcf.md (Status bump + weight-channel subsection),
forest-combiner.md (a stale "flooring" description corrected to the
snap), model-space-survey.md door 2 (shape (2) SHIPPED; the
non-interchangeability fact recorded), TODO, c-api-growth.md (the
dbarts_sampler_setForestWeights reservation with its three constraints),
NEWS. Sum-to-one unification landed on all three sites (verified live at
R_interface_bartcore.cpp:1158/:1330 pre-edit, R/A_class.R:402; the
shared constant sumToOneTolerance = 0x1p-26 lands at :63) with band-edge
tests per surface (1e-9 accepted, 1e-7 refused) in
test-sum-to-one-tolerance.R, isolating each site by direct slot mutation
past unrelated R-side gates. Gates double-run (implementer +
orchestrator): install, tinytest 3787/0, trio bitwise (27/27 + 3x5 + 5x7
vs c820227), air 0, lintr 0 on all three touched R files. Accepted
deviations: a create() helper to keep air-formatted .Call boilerplate
readable; the R validity test adjusted to dbartsModel's actual
initialize(proposal.probs=) signature.

ARC COMPLETE 2026-08-10: S0 4979656, S1 c820227 + d53efe6 (baseline
bcf-equivalence-c820227.rds), S2 153d1dd, S3 8c04f6e. All six CI
workflows green through the S2 engine push at landing time. Doors that
remain open are listed in "Doors held open"; the sum-to-one door is
CLOSED (unified at S3), and monotoneTreeIsFeasible's constant stays by
the recorded recommendation.

The reserved flat entry, dbarts_sampler_setForestWeights (S3's
docs/plans/c-api-growth.md reservation above), LANDED at dbarts-h-reshape
S1, ab3aa2fa, 2026-08-13; see that plan's S1 landing note and
c-api-growth.md's now-BUILT entry, not restated here.
