# multiforest-mutation-gaps

status: IN PROGRESS (commit 1 of 4 LANDED d489986 2026-08-06; items 2/3
        verified ready, H2 ready with its blocking clause bound; see the
        queue below)
agent: opus design + independent refuting critique with counterfactual
       power (two verification rounds), serialized opus implementers
rng: neutral throughout (refusals and opt-ins only; every landing gated
     on the equivalence trio bitwise)
budget: ~110 engine/bridge lines + ~165 test lines across 4 commits

## Goal

Close the recorded remainders of the multi-forest mutation surface
(TODO entry; doors recorded in docs/plans/runsbcbcf-repair.md Survey
addendum and Design + critique): multinomial setOffset/setTreatment
silent no-ops, the vacuously guarded BCF setWeights, and the
setState/installForests leaf-scale sibling door - plus the holes the
design census and critiques surfaced on the way.

## Process record

Design memo (opus, working probes in scratch clones) -> blind refuting
critique (opus, rebuilt every probe from scratch) -> focused
re-verification of the revised designs. The critique refuted item 3 as
designed (the (fitMin, fitMax) predicate is orthogonal to the leaf
calibration it protects: scaledResponseSd runs AFTER rescale(), so
BCF's s is affine-invariant in y - the guard admitted the divergent
same-endpoints/different-sd case and refused a harmless 10*y
cross-range restore), found the blocking H1 hole in item 1, refuted
the memo's D3 drift claim (installForests IS multi-forest-reachable:
test-interactions.R exercises BCF -> BCF), and found H2 out of scope.
The verification round then caught the H2 heteroscedastic blocker
(gaussian family + sigmaIsFixed_) and a WEAK-GATE hole in the item-3
test plan (tests/cpp statesAgree does not compare leafScale). Full
memo + critique texts were session-scratch (disposable); this file is
the durable record.

## Commit 1: offset/treatment refusals - LANDED d489986 2026-08-06

Multinomial setOffset and setTreatment were silent no-ops (verified
bitwise): the softmax is invariant to a common per-observation shift,
so a flat double* offset points along the null direction and a
meaningful one would be n x K - refusal is the correct surface, and
MultinomialResponse has no offset channel to begin with. Also closed
in the same commit, per the design census + critique:

- D1 (fourth hole, unrecorded): BCF setOffset(updateScale = TRUE)
  decalibrated the per-forest leaf scales exactly as the refused
  setResponse(TRUE) does (measured range 5.771 -> 21.176 with leaf
  scales unmoved). The offset conduit now carries the same two
  conditions as setResponse: multi-forest requires the combiner
  opt-in (supportsResponseMutation) AND updateScale exactly FALSE
  (NA refuses). setOffset(NULL, FALSE) still clears. The SUR-BCF /
  residual-conditioning embedding (correlated-outcomes.md) is
  untouched; the rbart_vi-style warmup-rescale convention is
  foreclosed for a BCF host, recorded there with the rejection of the
  recompute-s alternative (a data-dependent prior refresh mid-run).
- setTreatment now gates on the bcfGlue capability probe instead of
  numForests < 2, which a K-forest multinomial defeated; the
  single-forest message is preserved verbatim.
- D2: parseMultinomialData refused case weights and a test offset but
  silently dropped a host TRAIN offset; now refused at creation.
- H1 (critique, blocking): a multinomial TEST offset was refused at
  creation but accepted post-creation (bartcore_setTestOffset and
  setTestPredictorAndOffset gate only on refuseBCFTestSurface, which
  multinomial passes), and storeSample adds it to the K channels
  AFTER the softmax blend - reported probabilities left the simplex
  (measured per-row sums of 10). Both entries now refuse via
  refuseMultiForestTestOffset, ordered after refuseBCFTestSurface so
  BCF keeps its message; a null-offset test-predictor swap stays
  permitted.
- D4 + the H1 flat-API sibling: dbarts_sampler_setOffset and
  _setTestOffset carry the defensive refuseMultiForestMutation their
  setResponse/setWeights siblings got at 6744aca (unreachable today;
  symmetry).
- The stale MultinomialResponse comment (model.hpp) rewritten to
  cover both no-ops.

Gates, implementer + orchestrator re-run independently from fresh
--preclean installs: tests/cpp from clean all pass; tinytest
3616/3616 (+12); trio bitwise (27/27, BCF 5x6, multinomial 3x5);
sbc.R bcf 3 10 1 15/15 PASS; air clean; R CMD check Status: OK. No
dbarts.h change. Single-forest regression pins added, including
setOffset(TRUE) (rbart_vi's conduit).

## Queue (verified, not yet implemented)

2. Open BCF setWeights under the existing opt-in (design CONFIRMED by
   critique end to end): the guard is provably vacuous - swap-to-w2
   is bitwise create-with-w2 under resid.prior = fixed(); the
   default-prior sigest difference is identical to the shipped
   single-forest semantics (measured both). ~10 bridge lines reusing
   shape.supportsResponseMutation (no rename), comment extensions on
   the predicate, dbarts_sampler_setWeights stays refusing. Tests:
   flip the seam refusal, bitwise create-vs-swap arm under fixed(),
   fit.scale non-movement pin, zero-weight swap (numPositiveWeights_
   recount), multinomial still refused, single-forest untouched.
3. leaf.scale rides the state (INSTALL, the revised design, verified
   READY): append-only optional per-forest block (registry rule
   permits it, no stateFormatVersion bump), installed at BOTH
   setState and installForests - the state already restores k, the
   other half of the same leaf prior, and installForest already
   adopts the donor's transform/sigma/k/DART/glue. Verified: the
   divergent same-endpoints/different-sd case closes bitwise; the
   cross-range restore stays admitted; old states without the block
   reproduce today's behavior; downgrade (stock build reads a
   leaf.scale state) and upgrade both accepted through RDS; the
   single-forest change (donor node.scale now survives restore like k
   already does) completes an existing asymmetry and is correct.
   Consequence to record in the landing: a post-storeState
   setModel(node.scale) no longer survives save/load re-creation.
   Amendments bound: tests/cpp statesAgree must compare leafScale
   (WEAK-GATE otherwise); hostile-value posture matches k (non-finite
   /nonpositive read as absent or flow through; do not claim refusal).
4. H2, setSigma on models whose sigma is not free (found by the
   critique; higher severity than the three items - reachable through
   the SUPPORTED R5 sampler$setSigma): Chain::setSigma writes sigma_
   unconditionally and the constant-leaf draws divide by sigma_^2,
   but probit/logistic/multinomial pin sigma_ = 1 by the model
   definition (measured divergence 0.90-5.34) and the heteroscedastic
   variance forest owns the noise scale (measured 1.591). Guard at
   the bridge beside refuseBinaryWeightChange with the predicate
   family != gaussian && family != aft || shape.hasVarianceForest
   (the heteroscedastic clause is the verification round's BLOCKING
   amendment: buildVarianceForest sets sigmaIsFixed_ = true with
   family still gaussian). sigmaIsFixed_ alone is WRONG: gaussian +
   resid.prior = fixed() is the documented outer-Gibbs conduit
   (stan4bart mvbart drives it; breaking-with-lockstep was offered
   and declined - the family predicate wins on semantics, the
   gaussian-fixed path is intended behavior, not a grandfathered
   bug). Chain::setSigma itself must stay unguarded (setState
   restores through it). Includes the \item{sigma} .Rd correction
   (the length-1 broadcast, not one-per-chain).

## Related, recorded elsewhere

- The mu[leafOf]-gather SIMD item (docs/plans/setpredictor-leafof-
  rebuild.md door) ran its own design + critique in this session;
  record lands with that commit.
- Out of scope, left open: n x K multinomial offset (new shape, no
  demand); the setSigma-on-multinomial door subsumed into H2.
