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

## Status (2026-07-07)

Landed (933eed8). Chain::installForest seeds live trees, sigma, and k
from a donor's flat state - rng, latents, group effects, and saved
buffers stay fresh (a starting position, not a continuation) - via
rebuildLiveForest, extracted from setState so both share the rebuild
loop. Sampler::installForests cycles or user-maps donor samples onto
chains (NULL spreads them for overdispersed starts) and refuses
tree-count, cut-grid, and DART mismatches by name. R surface: the R5
installTrees(donor, samples = NULL) method plus bart2's warm.start
argument; the donor may be a sampler, a keepSampler fit, or a stored
state. Cross-grid remap REFUSED for v1 ("cross-grid warm starts are
not supported") - the stretch step stays open here. Deviations: the
diff ran ~630 insertions vs the ~250 budget - the bridge state parser
(mirroring setState's reader) and the test file account for it,
reviewed line-for-line as in-scope; step 3's component tests landed
as tinytest round-trips (getTrees/predict equality to 1.8e-14) since
getTrees is R-level; BCF glue transfers when both fits are BCF,
paralleling DART - a deliberate small superset of "sigma/k only".
Gates (implementer ran, reviewer re-ran): component tests pass,
tinytest 2469/0 (11 new), equivalence exact 18/18, air/lintr/checkRd
clean. Reviewer note (discharged 2026-07-07): rchk run over the 933eed8
tree, zero dbarts findings (7729 functions analyzed; only the usual
tool noise - SIMD flag warnings and R-internals abstraction errors).
The bridge SEXP parser is clean.
