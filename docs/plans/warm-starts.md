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
The bridge SEXP parser is clean. Post-landing fix (2026-07-08): CI's
valgrind leg caught installForests leaking its pool/sampleMap vectors
on every refusal path - the code freed the donor state before the
Rf_error longjmp but not the mapping vectors; they now live in a scope
that closes first (rchk could not see this: it checks PROTECT balance,
not C++ heap lifetimes).

Cross-grid landed (2026-07-22). The refusal is lifted: a donor fit on a
different cut grid (different n.cuts, useQuantiles, or predictor values)
now warm-starts by remapping. Sampler::installForests detects the grid
difference and routes each chain through Chain::installForestRemapped,
which builds every donor tree against the donor grid (recovering its
split indices) under a scoped ScopedCutGrid swap of the shared store's
cutPoints/numCuts - structural only, the live observation codes are
never touched - then reverts and applies Tree::mapOldCutPointsOntoNew
onto the live grid, collapsing starved splits, repartitioning, and
rebuilding fits, exactly the applyNewData body setData uses. The
interaction/column-mask containment pre-checks run under the same donor
grid (the remap only collapses, so a donor feasible pre-remap stays
feasible), keeping the refusal transactional. numTrees, DART, and
BCF-glue mismatches still refuse by name; a structural predictor
mismatch (categorical/continuous disagreement, or a malformed donor
grid) still refuses as gridMismatch, with a reworded message. Heap
lifetime on the error paths: the grid swap is a stack RAII guard
(ScopedCutGrid) whose destructor restores the live grid on every return,
including the buildFromFlat-failure early return, so no swapped store
ever escapes; the existing bridge free-before-longjmp scoping is
untouched. rng-neutral: the same-grid path still calls rebuildLiveForest
verbatim (installForest's donorCutPoints defaults to null), so a default
fit and a same-grid warm start stay byte-identical. Gates: component
tests pass (new testCrossGridWarmStart: a single-tree box donor on a
9-cut grid installs onto a 1-cut grid, collapsing 10 leaves to 2 with no
empty leaf and a finite run); full tinytest green (test-warm-start's
former cut-grid refusal is now a success assertion). Deviation from the
plan's "setData remap applied at install": the remap needed the donor
grid installed for the structural build, which only the sampler-owned
store can do - hence ScopedCutGrid and the mutable-store hand-off to the
chain, rather than the chain swapping its own const data_ reference.
