# warm-starts

agent: opus
rng: neutral engine-wise (a different, valid starting state; no
     default behavior changes)
budget: ~250 lines

## Goal

A sampler can start its chains from another fit's forest instead of
cold roots: an R-level verb installs saved trees (with leaf
parameters) as the initial state, for chain-count scaling, embedding
workflows, and - if grow-from-root lands - XBART-init-then-BART-refine.

## Context

- The machinery exists unexposed: flat trees replay through
  buildFromFlat, setState installs full chain states, and setData's
  Tree::mapOldCutPointsOntoNew remaps split indices across cut grids
  (src/bartcore/tree.hpp) - warm-starting onto a same-shaped sampler
  is state surgery, and onto a different cut grid is the setData
  remap applied at install.
- v1 scope: same predictors, same numTrees, donor fit has keepTrees
  or a stored state; each destination chain takes a donor sample
  (cycled or user-mapped). Cross-grid (new data) rides the existing
  remap machinery and is the stretch step.
- Everything downstream (fits, totalFits, sigma init) rebuilds through
  the existing state-install paths; no new invariants.

## Constraints

- Statistical honesty in docs: warm starts bias early samples toward
  the donor - burn-in guidance stays (shorter, not zero); say so in
  the Rd.
- numTrees mismatch refuses (subsetting/padding forests is a
  modeling decision, not a convenience).
- Out of scope: warm-starting hyperparameters beyond sigma/k (DART s
  transfers only if both fits use DART - refuse mismatches);
  cross-package tree imports.

## Steps

1. Engine: installForest(chain, flat trees + params) - a trimmed
   setState that touches trees/fits only, leaving RNG and
   hyperparameter state fresh; cut-grid remap path when the
   destination grid differs.
2. Bridge + R5 method + a bart2 `warm.start = <fit or sampler>`
   convenience (name at review); Rd with the two workflows.
3. Component tests: install-then-getTrees round trip; install from a
   fit onto identical data continues with well-formed trees and
   matching totalFits; cross-grid install collapses starved splits
   like setData does.
4. tinytest: end-to-end warm-started bart2 converges (statistical
   check) and beats cold start on early-iteration train RMSE for a
   deep-tree config (the measurable claim).

## Verification

- Component tests; full tinytest; equivalence exact (defaults
  untouched).
