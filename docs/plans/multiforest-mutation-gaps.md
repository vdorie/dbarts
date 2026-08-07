# multiforest-mutation-gaps

status: LANDED (all four commits, 2026-08-06: d489986 offset/treatment
        guards, ae7da51 setWeights opt-in, a7eb03d leaf.scale state
        block, ef61b39 setSigma guard; merged-tree battery green - see
        the landing notes)
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

## Commit 2: setWeights opt-in - LANDED ae7da51 2026-08-06

The multi-forest setWeights refusal was provably vacuous for BCF
(critique-confirmed end to end, both probes rebuilt independently):
Chain::setWeights is a pointer swap plus the numPositiveWeights_
recount, the combiner re-derives every per-forest response and
precision from y and w each sweep, and scaledResponseSd is unweighted,
so under resid.prior = fixed() swap-to-w2 is bitwise create-with-w2;
the default-prior sigest difference is identical to the shipped
single-forest semantics (measured both). The weight conduit now rides
shape.supportsResponseMutation (no rename; predicate comments on the
combiner/chain/shape now name the whole response-side conduit);
dbarts_sampler_setWeights stays refusing (flat-API precedent). Tests:
the seam refusal flipped, bitwise create-vs-swap arm under fixed()
with expect_identical on train/varcount/glue/both forests' fits plus
a not-inert sanity arm, fit.scale non-movement pin, zero-weight swap,
multinomial still refused, single-forest pinned. +112/-26; tinytest
3628 on its base.

## Commit 3: leaf.scale rides the state - LANDED a7eb03d 2026-08-06

INSTALL design (the critique refuted the original (fitMin, fitMax)
guard: scaledResponseSd runs after rescale(), so BCF's calibration is
affine-invariant in y and orthogonal to the transform - the guard
admitted the divergent same-endpoints/different-sd case, 0.335 on
fits of 0.82, and refused a harmless 10*y restore that agrees to
5.3e-15). The scale is an append-only optional per-forest block,
leaf.scale, after k in the registry (no stateFormatVersion bump),
written for every forest, installed by setState AND installForest
when positive; absent/non-positive/non-finite leaves construction's
scale (k's posture), so old states restore exactly as before, both
directions through RDS. The implementer found BOTH prior probe builds
had left the installForest arm dead: Sampler::installForests hands
installForest a reassembled ForestStateData, so the field must be
copied there (sampler.hpp) - the covering warm-start test was
falsified (commenting the copy fails it). statesAgree now compares
the field (the WEAK-GATE amendment). Consequences recorded in code
comments: a post-storeState setModel(node.scale) no longer survives
save/load re-creation (k's existing wart, applied to both halves of
the leaf prior); single-forest donors now carry node.scale through
setState as k already carried. The divergence closure is pinned
bitwise in test-bcf.R with a non-vacuity arm showing the refuted
transform guard would have admitted the case. +388/-3; tinytest 3639
on its base; ASAN tests/cpp clean.

## Commit 4: setSigma guard - LANDED ef61b39 2026-08-06

H2 (found by the critique; reachable through the SUPPORTED R5
sampler$setSigma): Chain::setSigma installs unconditionally and the
constant-leaf draws divide by sigma^2, but probit/logistic/
multinomial/ordinal/nbinom pin sigma at 1 by the model definition
(measured divergences 0.90-5.34) and the heteroscedastic variance
forest owns the residual scale (1.591) - the verification round's
BLOCKING clause, since buildVarianceForest leaves family_ gaussian
and would slip a family-only predicate. refusePinnedSigmaChange
(bartcore_bridge, external linkage) refuses on hasVarianceForest or
family not in {gaussian, aft}, with a separate accurate message per
branch; NEVER on sigmaIsFixed_ - gaussian + resid.prior = fixed() is
the documented outer-Gibbs conduit (stan4bart mvbart), pinned
positively by test. Chain::setSigma stays unguarded (both restore
paths reinstall the donor's sigma through it; commented in place).
The flat C API entry carries the same guard and it is LOAD-BEARING
there, not defensive: probit/logistic/ordinal/nbinom are flat-
creatable by family name and dbartsSpec(variance =) is an exported
consumer surface that reaches createHolder's variance control
(multinomial alone has no flat creation path). The \item{sigma} Rd
text now describes the length-1 broadcast and the refusal. All 16
in-tree setSigma callers verified gaussian; no assertion flipped.
+154/-2; tinytest 3624 on its base.

## Merged-tree battery (orchestrator, independent, at ef61b39)

install --preclean; tests/cpp from clean all pass; ASAN/UBSAN
tests/cpp leg all pass; tinytest 3659/3659 (exactly 3616 + 12 + 23 +
8); equivalence trio bitwise (27/27, BCF 5x6, multinomial 3x5); air
clean; R CMD check from a clean-copy tarball Status: OK. Commits 2-4
were each additionally gated by their implementers on their own bases
with the same battery (commit 2 also independently re-run by the
orchestrator pre-merge). sbc.R bcf smoke ran at commit 1; commits
2-4 do not touch the sbc paths (verified by reading in the critique's
gate check).

## Related, recorded elsewhere

- The mu[leafOf]-gather SIMD item (docs/plans/setpredictor-leafof-
  rebuild.md door) ran its own design + critique in this session;
  record lands with that commit.
- Out of scope, left open: n x K multinomial offset (new shape, no
  demand); the setSigma-on-multinomial door subsumed into H2; BCF
  whole-data setData, surveyed and held OPEN 2026-08-07 at the
  joint x/y/z-with-n-free semantics - record in
  docs/plans/runsbcbcf-repair.md ("setData door survey").
