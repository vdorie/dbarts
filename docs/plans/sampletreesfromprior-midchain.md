# sampletreesfromprior-midchain

agent: opus (design memo, refuting critique, implementer - serialized)
rng: neutral for every current call site in the repo and all four
  consumer repos; draws change only for mid-chain callers, which
  shipped code does not contain. Gate: the full equivalence trio
  bitwise, no snapshot moves.

## Defect

Found 2026-08-07 by the fused-suffstat design pass (recorded 6a94ec2).
Mid-chain sampleTreesFromPrior (chain.hpp:1203) zeroed muByTree to the
new arena while totalFits kept the pre-reset fits, so every later
sweep's residual rolled displaced - permanently, since
finalizeTotalFits overwrites from a residual that carries the term
(measured frozen to 8.9e-16 across six sweeps). On an arena shrink,
leafOf kept old node ids and rollTreeResidual read mu past .size()
inside capacity (ASAN-invisible under detect_container_overflow=0).
The invariant died at 917a53c, which added the monotone-required
muByTree resize without moving the map or the totals. Exposures: the
R method, the shipped ABI entry dbarts_sampler_sampleTreesFromPrior,
the setPredictor/setCutPoints transaction (revalidateTrees recovers
the zeroed parameters, so the pre-reset total survives; setData and
forceRefreshTrees self-heal), and growFromRoot after a reset (rolls
before rebuilding, never consults leafOfStale; measured 0.43
displacement vs 8e-16 control). Every existing call site was safe by
ordering luck: init-time buffers are already zero, and every repeated
site (samplePriorPredictive, the SBC harness) pairs the reset with
sampleNodeParametersFromPrior, which rebuilds all four buffers before
any read.

## Shape

Shape A narrow: the reset lands the forest in the ZERO-FIT state a
freshly built chain carries - totalFits/totalTestFits zeroed per
forest, vector/function treeFits rows zeroed, the constant-leaf map
memset to all-root (in bounds for every arena, gathers mu[0] == 0) -
with two deliberate deviations from a fresh chain, both load-bearing:
leafOfStale stays 1 (clearing it makes the fused suffstat pass
eligible one sweep early and would move sweep 1 of every default
bart2 run, whose init is this entry with no parameter draw after) and
muByTree is sized to the drawn tree (monotone reads it during the
first move). The reset is FOREST-ONLY: latents, sigma, k, DART, the
variance forest and the BCF/multinomial glue are untouched; a true
restart on a latent family follows with setResponse. Shape B
(init-only guard) is dead: every candidate predicate refuses the SBC
harness, whose rep-r+1 reset immediately follows rep r's run.
Design artifacts (memo, refuting critique, synthesis, consumer
census) are durable at <repo>/.claude/stfp-midchain-design/;
the critique verdict was STANDS WITH AMENDMENTS, its A1-A9 adopted -
notably A3 (the growFromRoot exposure the memo had declared safe) and
A1 (the displacement is totalOld minus a per-observation mixture of
stale OOB reads, not a clean carry; regression R^2 0.632 against the
clean-carry prediction of 1).

## Landing

LANDED 1947b10 2026-08-08 (squash of wt/stfp-midchain, base 49b36ad).
chain.hpp +93/-14: the reset tail, the rewritten contract comment
(the old one's "run() tolerates" claim had been false since 917a53c),
NDEBUG bounds asserts on leaf and leafPrev in rollTreeResidual, and
the leafOfStaleForTesting hook. tests/cpp +102 (testPriorResetContract:
the stale-flag pin T0, zero-totals and INV-1-after-a-sweep T1, the
flattened-prior arena-shrink probe T2 - reads only leafOf and
nodes.size(), never mu). New inst/tinytest/test-sampler-prior-midchain.R
(84 lines): train-vs-test falsifiers for gaussian, growFromRoot-after-
reset, and probit arms, each with its no-reset ulp control. Falsifier
discipline held: every behavioral test failed on the unpatched engine
(gaussian 4.24, growFromRoot 4.32, probit 0.277 vs ~1e-15 controls);
T0 plus the pre-existing fused-decline check (test_moves.cpp:1021-1022,
the trap detector) failed on the rebuildLeafOf trapped variant.
Gates, implementer and orchestrator independently: tests/cpp plain and
ASAN/UBSAN from make clean all-pass; tinytest 3665/0 (3659 + 6), no
snapshot moved; equivalence 27/27 "identical draws", BCF 5x6 and
multinomial 3x5 bitwise, no max-|z| line anywhere; air clean.

## Out of scope, recorded

- The setState round trip clears leafOfStale (rebuildLiveForest ->
  rebuildLeafOf) while a live post-reset chain keeps it at 1: a
  restored copy is fused-eligible one sweep earlier than its source.
  Unreachable on the normal path (every reset+parameter-draw pairing
  clears the flag before capture); record only.
- Rd wording for sampleTreesFromPrior/growFromRoot (the forest-only
  contract) can ride any future doc pass; the engine comment carries
  the contract.
