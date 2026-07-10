# ppd-sigma-pairing

agent: sonnet
rng: neutral for the sampler - pure post-processing of stored draws, no
  sampling-path changes, so equivalence is byte-identical for every
  scenario. ppd draws themselves DO change for multi-chain fits when the
  caller asks for combined output (extract/predict/fitted with the default
  combineChains = TRUE, or a fit whose own stored default is combined) -
  same RNG stream, corrected sigma pairing. Single-chain gaussian ppd draws,
  the chains-split gaussian path, and every binary ppd path are unchanged.
budget: ~30 lines net in generics.R, ~150 lines of tests, one plan file.

## Goal

Fix the defect pointwise-loglik surfaced and filed: for a multi-chain fit,
sampleFromPPD (R/generics.R) recycles object$sigma directly against ev.
object$sigma is always stored chain-interleaved (chain-fastest - draw s's
chain 1 value, then its chain 2 value, then draw s+1's), regardless of which
of its own two storage layouts a fit uses. Combined ev, however, is
chain-blocked (all of chain 1's draws, then chain 2's). Plain recycling
pairs sigma draw k with ev row k, which is the wrong chain past the first
n.samples rows - the noise added to a combined multi-chain gaussian ppd draw
uses another chain's sigma. Chains-split ev and single-chain fits were
already correct (verified: their row order matches sigma's own).

## Design decision

Two fixes were on the table:

(a) Reorder object$sigma from interleaved to chain-blocked before recycling
    it against combined ev, inside sampleFromPPD, keying off n.chains.
(b) Restructure the four ppd call sites to always hand sampleFromPPD
    chains-split ev (mirroring pointwiseLogLikelihood), then reshape the
    ppd result to the caller's requested layout afterward.

(a) was tried first and is statistically correct (each combined-row's noise
does use its own draw's sigma), but it does NOT give bit-identical draws
between a combined and a split ppd call under the same seed: reordering the
*sigma* fed to rnorm still draws its underlying uniform/normal variates in
whatever order the output array's flattening puts them in, which differs
between a chain-blocked matrix and a chain x sample x obs array. Verified
empirically - only the very first element coincides.

Landed instead: a third option, (a'), which gets (b)'s bit-identical
cross-layout guarantee without touching any of the four call sites. Inside
sampleFromPPD, the noise is always drawn in a single canonical order -
object$sigma's own chain-fastest order, spanning all n.chains * n.samples
draws per observation - regardless of what shape ev is in. When ev is
combined (chain-blocked, multi-chain), the freshly drawn noise vector is
reshaped - not redrawn - into that same row order using the package's own
combineChains() (bart.R), the identical helper the stored draws themselves
go through when a fit's own storage is combined. Because combineChains() is
a pure, lossless reshape, `combineChains(ev.split + noise.split) ==
ev.combined + combineChains(noise.split)` by linearity, so a same-seed
combined and split ppd call agree bit-for-bit after accounting for row
order - the same guarantee (b) would have given, at a fraction of the
diff, with no change to any caller.

n.chains is threaded into sampleFromPPD (a new 4th argument, default 1L)
since it is needed to build the reshape's target array dimensions; all four
callers (predict.bart, extract.bart, predict.rbart, extract.rbart) already
compute n.chains locally, so this is a one-line addition per call site.

## Implementation

R/generics.R, sampleFromPPD only:

- Both non-binary (gaussian) branches - unweighted and weighted - now draw
  a flat `noise` vector of length n.obs * n.draws in object$sigma's native
  order (unweighted: `rnorm(n.obs * n.draws, 0, rep_len(object$sigma, ...))`;
  weighted: same sd vector as before, `sigma / sqrt(w)` per observation
  block, order otherwise unchanged).
- When `n.chains > 1L && length(dim(ev)) < 3L` (ev is combined), the noise
  vector is reshaped via `combineChains(array(noise, c(n.chains, n.draws %/%
  n.chains, n.obs)))` before being added to ev. Otherwise (chains-split or
  single-chain) it is added as the same flat vector the pre-existing code
  used - bit-identical to before.
- Binary branches (rbinom, weighted and unweighted) are untouched: they use
  no sigma, so there was never a pairing defect to fix, and this stays true
  regardless of ev's shape.
- object$seed save/restore at the top and bottom is untouched.

Doc comment above sampleFromPPD rewritten (the old one was already flagged
"outdated" by the pointwise-loglik landing) to state the new invariant.

Call sites: predict.bart, extract.bart, predict.rbart, extract.rbart each
gain `n.chains` as the sampleFromPPD call's 4th argument - all four already
compute n.chains earlier in the same function for their own combineChains
handling, so no new lookups.

## Test

New file inst/tinytest/test-ppd-sigma-pairing.R (a dedicated home rather
than extending test-generics-posteriorPredictiveDistribution.R, which is a
general-purpose ppd smoke file predating this fix and stays that way):

1. Layout invariance: a 3-chain bart2 gaussian fit, same seed reset before
   two extract(type = "ppd") calls (combineChains = FALSE and the default
   TRUE); `dbarts:::combineChains(ppd.split)` must equal `ppd.comb` exactly
   (expect_identical), proving the cross-layout bit-agreement the design
   relies on.
2. Within-draw coupling: a short (12 samples, 5 trees), 4-chain fit on 15
   observations, chosen so sigma varies noticeably draw to draw. The
   combined ppd draw is recomputed by hand - rnorm in sigma's native order,
   reshaped via the same combineChains() call sampleFromPPD makes - and
   compared bit-for-bit (expect_identical). This is the regression test:
   verified to fail (checked out the pre-fix generics.R, reinstalled, reran)
   for every case but the untouched single-chain one.
3. Single chain: ev + rnorm(sigma) with no reshape, matching both the new
   and the pre-fix code exactly - confirms the fix leaves this path alone.
4. Weighted gaussian, 3-chain combined: manual recomputation includes
   sigma / sqrt(w) per observation on top of the same reshape, confirming
   weight recycling and the sigma fix compose correctly.
5. rbart_vi: a 2-chain fit, same layout-invariance check as (1) through
   extract(type = "ppd").

inst/tinytest/test-generics-posteriorPredictiveDistribution.R had a
pre-existing multi-chain combined check (added before this fix existed)
that encoded the old, incorrect pairing directly: `rnorm(n.samples *
n.chains, 0, bartFit$sigma)` recycled without any reorder or reshape. It is
a formula, not a hardcoded numeric literal, so it was updated in place to
recompute the same way sampleFromPPD now does (draw in native order, reshape
via dbarts:::combineChains) rather than left to encode the bug; the
chains-split block two paragraphs below it needed no change (already
correct, confirmed by rerunning against the pre-fix code). No other test
file hardcodes multi-chain combined gaussian ppd values.

## Verification

Gates, worktree at .claude/worktrees/ppd-sigma-pairing, private library
/Users/vdorie/.claude/jobs/7fe13675/tmp/Rlib-land:

1. R CMD INSTALL -l <lib> (R-only, no --preclean) - DONE (dbarts).
2. tinytest::test_package("dbarts") - 2590 passed / 0 failed (2581 baseline
   + 9 new checks).
3. Equivalence vs benchmarks/baselines/equivalence-de67cbb.rds - all 21
   scenarios "identical draws (same RNG stream)", 21 compared / 0 skipped.
4. air format --check and lintr on all three touched R files
   (R/generics.R, inst/tinytest/test-generics-posteriorPredictiveDistribution.R,
   inst/tinytest/test-ppd-sigma-pairing.R) - clean, 0 lints.
5. Pre-fix regression check: reverted generics.R to HEAD (git stash),
   reinstalled, reran test-ppd-sigma-pairing.R - checks 1, 2, 4, 5 (every
   multi-chain case) failed as expected; check 3 (single chain) still
   passed, confirming that path is genuinely untouched. Restored the fix
   (git stash pop) before landing.

## Status

- 2026-07-10: implemented on wt/ppd-sigma-pairing; all gates green
  (tinytest 2590/0, equivalence 21/21 identical, 0 lints on all three
  touched files, pre-fix regression check confirmed the new tests catch
  the bug in every multi-chain case and leave single-chain alone).

- 2026-07-10: LANDED as 03a5b85 (squash of wt/ppd-sigma-pairing).
  Reviewer verified the pairing arithmetic in all four branches
  (native-order noise + combineChains reshape) and re-ran gates on
  the landed shared-library build: tinytest 2590/0, equivalence
  21/21 identical. The implementer's pre-fix regression run (4/5
  scenarios fail on the old code, single-chain passes) is the
  defect's closing proof.
