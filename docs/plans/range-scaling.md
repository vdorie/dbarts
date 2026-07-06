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
