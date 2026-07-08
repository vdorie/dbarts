# moves-degenerate-root-guard

agent: opus
rng: neutral (guard fires only on the degenerate case that crashes today;
     non-degenerate trees keep their exact code path and RNG stream)
budget: guard ~10 lines; tinytest + tests/cpp cases ~70 lines

## Goal

A root-only tree with no available split variable stops crashing the
sampler. Birth/death treats such a tree's move as a no-op for the sweep
(no birth, no death; the leaf draw still happens as usual), instead of
forcing a birth that draws a rule for invalidVariable and indexes
data.types out of bounds.

## Context

- Defect: correctness-audit block-1 finding 1. drawBirthableNode's
  hasSingleNode branch (moves.hpp ~L104) returns the root with prob 1.0
  without checking availability; probabilityOfBirthStep (~L64) then
  returns 1.0 for the single node, so birthOrDeathMove (~L138) forces a
  birth and drawRuleAndVariable -> drawSplitVariable returns
  invalidVariable (-1) -> drawRuleForVariable reads
  data.types[(size_t)-1] (model.hpp ~L1583). The death branch is
  equally unguarded (fillNoGrand is empty for a single node).
- Latent on normal data (a constant column still quantizes to >= 1 cut;
  those births hit the empty-leaf veto and reject legally), reachable
  via zero-cut columns: setCutPoints with empty vectors on a root-only
  sampler (sampler.hpp setCutPoints; invalidCutPoints only guards
  existing splits).

## Constraints

- RNG-neutral: non-degenerate data must draw the identical stream
  (equivalence exact 18/18).
- Engine house style: ASCII, LLVM-Doxygen comments only where the code
  cannot speak; comment delta <= code delta.

## Steps

1. Guard the top of birthOrDeathMove: if hasSingleNode() and no
   birthable node exists, set stepTaken/stepWasBirth false and return.
   birthableNodeExists is only reached for single-node trees, and a
   movable tree never triggers it, so no live stream shifts.
2. tinytest: setCutPoints to zero-length vectors on a root-only
   sampler, run, assert finite output and root-only trees.
3. tests/cpp: construct the degenerate state directly (private rng +
   deterministic data so the shared streams do not shift), run, assert
   sigma finite and trees still root-only.

## Verification

- Reproduce the crash pre-fix (separate Rscript subprocess), then the
  same script clean post-fix; tests/cpp; full tinytest; equivalence
  exact 18/18.

## Status (2026-07-08, LANDED)

Reproduced pre-fix: the unfixed package segfaults at sampler$run
("caught segfault ... invalid permissions", exit 139) on a root-only
sampler after setCutPoints(list(numeric(0), numeric(0))). Post-fix the
same script runs clean: predictions finite, max nodes per tree 1 (all
leaves), sigma finite, exit 0.

Guard shape: at the top of birthOrDeathMove,
`if (tree.hasSingleNode() && !birthableNodeExists(ctx, tree)) {
stepTaken = stepWasBirth = false; return 0.0; }`. Short-circuits on
hasSingleNode so multi-node trees never call birthableNodeExists;
single-node trees that can split evaluate it to true and proceed
unchanged - no RNG consumed either way.

Gates: R CMD INSTALL --preclean clean; tests/cpp all tests passed (new
"degenerate root guard" case); tinytest 2473 ok / 0 fail; equivalence
exact 18/18 identical draws. The C++ case uses its own rng and
deterministic data so the shared ext_rng and global runif01 streams,
and thus the downstream model/gp suites, are unperturbed.
