# empty-leaf-veto

agent: opus
rng: posterior-changing (proposal kernel changes; stationary
     distribution preserved but gated as the strictest class)
budget: ~200 lines + a short design note

## Goal

The no-empty-leaf invariant is enforced one way, not two: ordinal
proposals become occupancy-aware (as categorical rules already are by
gauge construction), or - if the investigation says the cost is not
worth it - the -1e7 veto is kept but documented as a deliberate,
finite, load-bearing constant.

## Context

- Veto: src/bartcore/moves.hpp:57-59, ported verbatim from classic
  likelihood.cpp. Categorical rules cannot produce empty sides
  (canonical gauge, docs/design/core-generalization.md); ordinal births
  draw occupancy-blind and rely on the veto.
- The finite value is load-bearing: vetoed-vs-vetoed branch comparisons
  (both sides -1e7) reduce to the prior ratio, which is exactly how the
  gp-leaves over-cap deadlock was diagnosed (docs/design/gp-leaves.md,
  superseded veto precedent). -inf would produce NaN there. Any change
  must preserve defined vetoed-vs-vetoed behavior or make it
  unreachable.
- Occupancy-aware birth: a node's valid ordinal cut set is
  [min code, max code) over its segment (left = code <= cut); min/max
  is one O(n_leaf) scan. Restricting the proposal changes its density,
  so the birth/death MH ratio needs the matching reverse-proposal term.
  change/swap validate descendants via the good-rule flows but still
  rely on the veto for occupancy - full removal must cover them too.

## Constraints

- Same stationary distribution: the exact-posterior gates are the
  arbiter, not opinion.
- If the full occupancy-aware move set (birth + change + swap) exceeds
  budget, stop after the investigation step and write the design note
  recommending keep-and-document; that is a valid outcome.
- Out of scope: gp-leaves' over-cap fallback (separate mechanism).

## Steps

1. Investigate: enumerate every read of the -1e7 constant and every
   move path that can create an empty leaf; write up whether
   vetoed-vs-vetoed is reachable outside gp over-cap.
2. If proceeding: occupancy-aware ordinal birth (min/max scan, adjusted
   proposal densities both directions); extend the change-move
   rejection sampler and swap validity walk to occupancy; delete the
   veto.
3. Design note in docs/design/ recording whichever outcome.
4. Regenerate snapshots; re-record equivalence baseline.

## Verification

- benchmarks/R/logistic-reference.R and categorical-exact.R match the
  exact posterior to MC error.
- Component tests (move reversibility on fixed inputs); full tinytest.
- Equivalence z-mode vs prior baseline; bench-sampler compare (the
  min/max scan adds per-birth work).
