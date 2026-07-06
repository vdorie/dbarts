# range-scaling

agent: fable (decision memo); implementation only if overturned
rng: posterior-changing if changed
window: pre-release (changing the default scaling after release
        re-litigates every user's defaults)
budget: memo only; ~1 page added to docs/design/

## Goal

The response scaling (y mapped to [-0.5, 0.5] by range) is a recorded
decision instead of an unexamined BayesTree inheritance. The likely
outcome is keep-and-document; the memo makes that an argument rather
than an accident.

## Decision

Question: range scaling vs sd standardization for continuous responses.

For keeping range scaling: it is the convention of the entire BART
software lineage (BayesTree, BART, bartMachine), so priors and k values
transfer across packages and papers; the k = 2 / node.scale = 0.5
defaults are calibrated against it; changing it invalidates every
cross-package comparison and all RNG-locked artifacts at once.

Against: outlier sensitivity (two extreme y values compress the
effective prior on everything else); it is the root cause of the
internal-scale bookkeeping (sigma stored internal-scale,
sigmaPriorScale, the restoreScale round-trips -
src/bartcore/model.hpp:1841-1889, chain.hpp:203-206), though
state-continuation removes most of that cost independently.

Recommendation: keep; document the choice and the outlier caveat in
man/dbarts.Rd and a paragraph in docs/design/prior-defaults.md
(see prior-constants). Evidence that would change it: a demonstrated
real-data failure mode where sd scaling materially outperforms, or a
decision to break with the lineage defaults generally.

Response (VD, 2026-07-06): before deciding, survey what other BART
implementations do; and the memo must address Gibbs-sampler drift -
rescaling on setResponse lets the effective prior drift during burn-in,
while proper samples need the scale locked to a fixed range.

Survey (2026-07-06, primary sources): range-anchored leaf priors are
the classic lineage (BayesTree, BART via tau = range/(2k sqrt(m)),
bartMachine, bartpy; flexBART z-scores y then range-anchors tau);
every post-2018 design anchors to sd/var(y) instead (XBART, stochtree,
bcf, SoftBart's newer interface, PyMC-BART), typically with a
hyperprior on the leaf scale as the robustness mitigation - which
dbarts already ships as the chi(1.25, Inf) prior on k. The only
published side-by-side is Linero (arXiv:2210.16375 sec 1.5): "I have
not found the choice of standardization for Y_i to have a large
impact on the results"; the only published outlier caveat is
bartMachine's one-liner (workaround: log/winsorize y). bcf is the
clearest reasoned departure, argument implicit.

Code finding: setResponse re-anchors the internal scale on EVERY call
(model.hpp:1764-1775 rescale(); classic did the same, bartFit.cpp
setResponse - bartcore inherited it faithfully), holding sigma and
the variance prior
fixed on the original scale but letting the leaf prior's
original-scale width track the new response range - VD's drift,
live today. setOffset re-anchors only with updateScale = TRUE
(documented burn-in-only); setResponse's behavior is undocumented.
No equivalence scenario covers mid-run setResponse.

Revised recommendation: (1) keep range anchoring (lineage k transfer;
k hyperprior is the robustness escape; document the outlier caveat);
(2) fix the semantics, not the anchor - setResponse gains updateScale
defaulting to FALSE (locked at creation), mirroring setOffset; Gibbs
users re-anchor explicitly during burn-in, then lock. rng: shifting
only for mid-run setResponse users; standard fits unchanged. setData
keeps always-rescale (structural replacement).

Signed off (VD, 2026-07-06): keep, document, add the setResponse
updateScale switch with default FALSE. bcf's sd anchoring is
considered only within forest-split-bcf. Documentation folds into
prior-constants; the switch is the remaining implementation:

## Steps (setResponse updateScale)

agent: opus; rng: shifting (mid-run setResponse tests only)

1. Engine: ResponseModel::setResponse gains updateScale; FALSE path
   mirrors setOffset's reuse-existing-scale branch (model.hpp:1814).
2. Bridge + R5 setResponse(y, updateScale = FALSE); dbarts.h stays
   (additive change only if a C consumer needs it - none known).
3. man/dbartsSampler-class.Rd documents both switches and the
   burn-in-only guidance; snapshot regeneration for affected test
   files; component test covering locked vs re-anchored.

## Steps

1. Memo above lands as the design-note paragraph; VD signs off on
   keep vs change.
2. If keep: prior-constants absorbs the documentation step; close.
3. If change: write a fresh implementation plan (re-derive
   node.scale/sigquant anchoring, full snapshot regeneration,
   exact-posterior gates, equivalence baseline break) - not this file.

## Verification

- The design note exists and man/dbarts.Rd states the scaling and its
  caveat. No code change under the recommended outcome.
