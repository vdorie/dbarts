# chi-k-empty-leaf-count

agent: opus
rng: neutral (the empty-leaf veto keeps MCMC moves from stranding a leaf,
     so every fit without a stranded empty leaf draws the identical stream;
     equivalence exact 18/18, setdata included)
budget: fix ~4 lines; tests/cpp case + chain test hook ~90 lines

## Goal

A leaf that a data mutation strands empty stops corrupting the chi-k
hyperprior's accounting. A forced-to-zero empty leaf is not a draw from
the k-scaled prior, so it must not enter the count or the sum of squares
the ChiKHyperprior draw consumes - matching the function-leaf path, which
already excludes it.

## Context

- Defect: correctness-audit block-5 finding (latent). In
  src/bartcore/chain.hpp sampleParametersAndSetFits, the constant-leaf
  branch gives an empty leaf param = 0.0 but still ran
  `kSumSquaredParams += 0` and `kNumLeaves += 1.0`; the vector (linear)
  branch kept its zero block but still ran `kNumLeaves += numParams`.
  Empty leaves thus inflated the posterior shape (0.5 (totalNumLeaves +
  df)) while contributing nothing to the sum of squares, biasing k.
- The function-leaf (GP) branch is correct already: model.hpp's
  drawFromPosteriorForNode / drawFromPriorForNode return
  FunctionLeafDrawStats{0, 0} for empty leaves, so the branch's
  `kSumSquaredParams += stats.sumSquaredParams; kNumLeaves +=
  stats.numParams` adds nothing.
- The ChiKHyperprior::draw itself (model.hpp ~L1727) is correct and
  unchanged. Its gate `kSumSquaredParams > 0.0` (chain.hpp ~L688) is
  unchanged.
- Latent because no public mutation strands an empty leaf: applyNewData
  (setData) and revalidateTrees (setPredictor, rejected on any empty
  leaf) both collapse empty nodes, and collapse bubbles up while n > 0,
  so a live sweep never sees one. Only fabrication reaches the branch.

## Fix

- Constant and vector branches: gate the k-stat accumulation on
  `node.numObservations() > 0`, so a forced-zero empty leaf is skipped in
  both count and sum of squares, matching the GP semantics. Grep confirms
  the only accumulation sites are these two plus the correct GP branch;
  no other site touches kNumLeaves / kSumSquaredParams.

## Verification

- tests/cpp regression testChiKEmptyLeafAccounting (test_model.cpp,
  beside testChiKHyperprior): a chain hook
  (Chain::accountStrandedLeafKStatsForTesting) births tree 0, strands its
  right child empty, and runs the real per-sweep parameter draw with
  updateK on, returning the k accounting. The test asserts the accrued
  leaf count equals the populated leaves only (stride 1 scalar, 2 linear).
  Fails on the unfixed code (counted the empty leaf: 2 vs 1, 4 vs 2);
  passes fixed. Both inline-accounting branches are exercised.
- The hook is the seam because the public mutation surface cannot strand
  an empty leaf (see Context).

## Status (2026-07-08, LANDED)

Fail-on-unfixed confirmed: with the two guards reverted, the rebuilt
tests/cpp reports `FAIL: constant leaf excludes the empty leaf` and
`FAIL: linear leaf excludes the empty leaf` (the fabricated
populated-plus-empty guard still passes); re-applying the guards clears
both.

Gates (from the worktree):
- R CMD INSTALL --preclean . : clean (DONE).
- tests/cpp: make clean && make && ./test_bartcore -> all tests passed,
  including "ok: chi-k empty-leaf accounting".
- tinytest::test_package("dbarts") -> 2465 pass / 0 fail.
- equivalence compare equivalence-cf99a00.rds -> 18 / 18 identical draws
  (same RNG stream), setdata included; "posteriors statistically
  indistinguishable". No divergence.
