# forest-cache-drift

Status: PROPOSED, 2026-09-24. Decision-gated: VD signs off D1 to D5 before implementation.
agent: opus (engine, component tests, gate battery, re-records); sonnet (records: docs, INDEX, TODO)
rng: posterior-changing. Every chain in which a forest travels the amplitude ridge draws differently
  from its first rescale on, and the implemented chain's stationary law moves, because the defect
  biased it. Every other configuration (single-forest, heteroscedastic, multinomial, and an amplitude
  coupling whose forests are all held or unridged) must replay its baseline bitwise.
window: none (VD 2026-09-24: the right fix, whatever its timing against the merge).
budget: engine ~40 lines (combiner.hpp, chain.hpp); tests/cpp ~200; tinytest ~10; benchmarks ~80
  (BCF bench scenarios, a pooled-seed gate mode); two baseline re-records; records.

## Goal

A forest's cached fits never carry a rounding gap that a transform has multiplied. Every forest that
travels the amplitude ridge ends each sweep with `totalFits` bitwise equal to its leaf tables gathered
in tree order. Every other cache differs from its derivation by accumulated additive rounding only, and a restore still
drops that rounding-level gap for forests the ridge does not move. The
latent BCF exact gate passes in both modes on both architectures, and the gate and calibration records
that read the drift as mixing are measured again.

## Context

The defect. [`AmplitudeForestCombiner::rescaleAmplitudeRidge`](../../src/bartcore/combiner.hpp) moves
forest f along its ridge, `(a_f, leaves) -> (a_f/c, c leaves)`, once per sweep, from
[`AmplitudeForestCombiner::afterCombine`](../../src/bartcore/combiner.hpp). That is reached from
[`Chain::run`](../../src/bartcore/chain.hpp), [`Chain::growForestFromRoot`](../../src/bartcore/chain.hpp)
and [`Chain::interweaveGlueRidgeForTesting`](../../src/bartcore/chain.hpp). It multiplies `muByTree` by
c, and separately multiplies the [`Forest`](../../src/bartcore/combiner.hpp)'s `totalFits` by c. Let
`g_i = totalFits_i - sum_t mu_t[leafOf_t(i)]`. The sweep keeps `totalFits` by difference updates
([`Chain::rollTreeResidual`](../../src/bartcore/chain.hpp),
[`Chain::fusedRollPass`](../../src/bartcore/chain.hpp),
[`Chain::finalizeTotalFits`](../../src/bartcore/chain.hpp)), which carry g unchanged through every leaf
redraw and add a rounding increment. The rescale maps g to `c g`. `a_f` is redrawn from its conditional
between rescales, so the c factors are not ratios of one quantity and do not cancel: g follows a
multiplicative random walk while the fits stay O(1).

This has two effects. First, g enters every tree's partial residual (tree t fits `y - g - sum_{s != t}
f_s`), so the leaves absorb `-g` while `totalFits` still tracks the data. The sigma draw (through the combined
location), the latent refresh, the amplitude conditional and the per-forest reporting channel then read one function (the
cache). The leaf prior, the ridge's GIG statistic M, the saved trees, predict and the state read another
(the leaves). Second, [`Chain::setState`](../../src/bartcore/chain.hpp) and
[`Chain::installForest`](../../src/bartcore/chain.hpp) re-derive the cache
([`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp)). So a restored chain drops g without warning,
and continues from a different state than the one that was saved.

Measured by the investigation's probes, which were not kept (Steps 1 and 7 replace them). At the exact
gate's single-tree shape, the index-scale gap `|a| max|fit - leaf|` grew from 1e-11 to 0.3 under probit
(seed 13). Within 20000 sweeps it reached 0.53 under logistic and 2.4e-3 under gaussian. Latent BCF tau
was biased +0.005 at the gate's quick length and +0.027 at full length, a pooled z of 7.1 over 303 seeds.
That bias is what [`exactLatentBCF`](../../benchmarks/R/bcf-latent-exact.R) caught in quick mode on arm64,
and it does not depend on the host. A scratch patch that rebuilt `totalFits` from the rescaled leaves
zeroed the gap and restored the gate.

Which models are affected: every forest that the bridge gives a half-Cauchy amplitude.
[`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp) sets `ridge` from `amplitudePriorScale > 0`.
That covers bcf's prognostic forest and, in the formula route's K-forest fits, each forest without a
basis (`forestParams` in [R/model.R](../../R/model.R); in bart2twoforest, forest 0 only), under gaussian,
probit and logistic. bcf's treatment forest ships with the ridge off
([ridgeB is code that is OFF](../design/multiplier-combiner.md#ridgeb-is-code-that-is-off)).

Design: [The ASIS ridge](../design/multiplier-combiner.md#the-asis-ridge), whose "rescale-consistency
set" lists `totalFits` as a quantity the move scales. This item reverses that.
[5. Correctness](../design/level-fibre.md#5-correctness),
[Identification and the level-centering move](../design/multinomial.md#identification-and-the-level-centering-move),
[Decision 2 - the exact gate](bcf-latent-evidence.md#decision-2---the-exact-gate).

## Audit

This covers every site in src/bartcore and the bridge that transforms leaf values, a calibration, a
response scale or a forest's level in bulk, and what each does to derived caches. "Amplifies" means the
gap is multiplied. "Preserves" means the gap is carried with additive rounding. "Clears" means the cache
is re-derived from scratch.

| site | transform | caches, as adjusted | gap | fix |
|---|---|---|---|---|
| [`AmplitudeForestCombiner::rescaleAmplitudeRidge`](../../src/bartcore/combiner.hpp) | leaves times c | `totalFits` times c, separately | amplifies, never cleared | YES |
| same, `record` branch | as above | `totalTestFits`, `currTestFits` times c | preserves until [`Chain::run`](../../src/bartcore/chain.hpp) rebuilds them from zero at the next recorded sweep; BCF reports no test channel | no |
| same, keepTrees slot | saved leaf values times c | the same multiply as the live leaves, so bitwise the live values | none | no |
| [`MultinomialForestCombiner::afterCombine`](../../src/bartcore/combiner.hpp) | occupied leaves plus c/m | `totalFits` plus c (m fl(c/m) is not exactly c); test fits untouched by design | preserves: additive, rounding-level | no (D2) |
| [`Chain::drawLevelShift`](../../src/bartcore/chain.hpp) | occupied leaves plus c_t, where the c_t sum to 0 | `totalFits` untouched, correct algebraically | preserves: the projection's rounding residual | no (D2) |
| [`Chain::reanchorVarianceForest`](../../src/bartcore/chain.hpp) | variance factors times g | `combinedVariance` and its test twin times f, where g^m' is not exactly f | preserves for one sweep; [`Chain::sweepVarianceForest`](../../src/bartcore/chain.hpp)'s closing product and [`Chain::refreshVarianceTestFits`](../../src/bartcore/chain.hpp) clear it | no |
| [`GaussianResponse::setResponse`](../../src/bartcore/model.hpp), [`GaussianResponse::setOffset`](../../src/bartcore/model.hpp), [`GaussianResponse::setData`](../../src/bartcore/model.hpp), [`GaussianResponse::restoreScale`](../../src/bartcore/model.hpp) | the response transform moves; leaves are reinterpreted, not transformed | none | none | no |
| [`Chain::setForestPriorScale`](../../src/bartcore/chain.hpp), [`Chain::setModel`](../../src/bartcore/chain.hpp), [`Chain::setForestBasis`](../../src/bartcore/chain.hpp) | leaf scale, k | none: no drawn leaf is touched | none | no |
| [`Chain::rebuildFitsFromParameters`](../../src/bartcore/chain.hpp), and the per-observation update session, which ends in the same rebuild | per-tree re-route, leaves unchanged | `totalFits` subtract, then add | preserves | no |
| sweep upkeep: the three roll functions above, [`VarianceForest::applyLeafFactor`](../../src/bartcore/chain.hpp) | per-tree redraws | difference updates | preserves (1.2e-15 at 200 trees, [`runEnsembleTests`](../../tests/cpp/test_ensemble.cpp)) | no |
| [`Chain::sampleTreesFromPrior`](../../src/bartcore/chain.hpp), [`Chain::sampleNodeParametersFromPrior`](../../src/bartcore/chain.hpp), [`Chain::applyNewData`](../../src/bartcore/chain.hpp), [`Chain::forceRefreshTrees`](../../src/bartcore/chain.hpp), [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp), [`Chain::rebuildLiveForestRemapped`](../../src/bartcore/chain.hpp), [`Chain::sampleVarianceForestFromPrior`](../../src/bartcore/chain.hpp), [`Chain::rebuildVarianceForest`](../../src/bartcore/chain.hpp), [`Chain::refreshVarianceForest`](../../src/bartcore/chain.hpp) | wholesale replacement | rebuilt from zero | clears | no |
| the bridge, [`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp) and the state entries | none: forest state is only passed through, and state goes through the re-deriving restore | reporting reads the cache (`$getForestFits`) and the leaves (trees, predict) as they are, and the gap made them disagree | none of its own | no |

The remaining caches are rebuilt every sweep or scratch, and no transform reaches them: node sufficient
statistics and `treeY` (rebuilt by the roll, dead across sweep boundaries), the multinomial offset slab
(rematerialized at sweep entry and at reporting), the linear leaf's U'WU cache (no bulk transform of vector
leaves exists), and combiner scratch.

Two constraints stand. First, `kSumSquaredParams` is accumulated from pre-rescale leaves and consumed
after `afterCombine`. That is unreachable only because
[`Chain::buildSpecifiedForest`](../../src/bartcore/chain.hpp) and
[`Chain::buildMultinomialForest`](../../src/bartcore/chain.hpp) pin `updateK` false. Opening a k
hyperprior on a ridged forest must move the accumulation after the rescale. Second, the TODO item
level-shift-capi would add the first transform a host drives. Its transform is additive, so rule (ii)
below covers it. Any host-facing multiplicative transform falls under rule (i).

## Decision

D1, the fix shape.

- (a) Re-derive the cache after the transform. One helper on `Forest` gathers `totalFits` from
  `muByTree` through `leafOf`, seeded at zero and summed in tree order. That is exactly the arithmetic of
  [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp) through
  [`Chain::addTreeFitsToTotal`](../../src/bartcore/chain.hpp), so a restore at a sweep boundary reproduces
  the cache bitwise. `rescaleAmplitudeRidge` stops writing `totalFits`. `afterCombine` calls the helper
  for every forest it moves (`update && ridge`), after the rescale and on its no-op returns too, so the
  invariant holds at every sweep boundary rather than only on sweeps where c moved.
  - Stale leaf maps are NOT refreshed. Only [`Chain::sampleTreesFromPrior`](../../src/bartcore/chain.hpp)
    marks a map stale. It gives the tree an all-zero leaf table over an all-root map, which holds until
    the tree's own draw rebuilds the map. The stale pair therefore gathers exactly the zero the cache
    holds for it. Refreshing would clear `leafOfStale`, which would make the fused pass eligible a sweep
    early and move draws; that is also why `drawLevelShift` declines stale trees rather than rebuilding
    them. Inside a sweep no map is stale by the time `afterCombine` runs anyway. The helper asserts the
    zero table under `!NDEBUG`. This replaces the scratch patch's fallback to the old scaling on a stale
    map, which kept the defect on exactly that branch.
  - Cost: per row per tree, a leaf-map read and a `totalFits` read-modify-write, about 20 bytes against
    the residual roll's 24, so close to the roll's own traffic. The pre-fusion gaussian x86 profile
    ([10. Re-profile and census at 06f73b0 (2026-08-04, dbarts-bench)](../design/memory-wall-frontier.md#10-re-profile-and-census-at-06f73b0-2026-08-04-dbarts-bench))
    put the roll at 25 to 28 percent of a sweep, and the fused pass has since made the sweep faster, so
    the roll's share is larger now. The critique timed both forms on a loaded arm64 machine (indicative
    only): the plain gather cost 7 to 13 percent of a BCF sweep, a blocked gather 6 to 11. The blocked
    form reads and writes `totalFits` once per block of trees while adding each row's trees in order, so
    it gives values identical to the plain gather (identical draws over 500 probit sweeps). The blocked
    form is the default; the block width is tuned in Step 6.
  - A check in debug builds only: at the end of every sweep of `Chain::run` under `!NDEBUG`, every
    moved forest must match the gather bitwise, and every other constant-leaf forest of a double-precision
    chain must satisfy `|totalFits_i - gather_i| <= 1e-8 (1 + |totalFits_i|)`. Float chains are exempt:
    the fp32 tests carry gaps above 1e-10 by design (`testEndToEndGaussianFp32`, and the fp32 arm of
    `testGatherTailShapes` checks at 1e-5). The logistic treatment forest, ridge off, reached 6e-11
    relative within 20000 sweeps with the fix in, so 1e-10 would have no margin. The R build
    defines NDEBUG and compiles the check out. It runs in tests/cpp, where it catches any later transform
    that breaks the rule below, in any test that reaches that transform.
- (b) A per-forest scale applied lazily. Not a fix: effective leaves are `s_f mu`, and the move sets `s_f *= c`,
  `a_f /= c` in O(1). Nothing stored is multiplied, so no gap is. It touches: `forestMultiplier`, and
  through it `formForestResponse`, `formForestVetoWeights`, `combinedFits` and `drawForestAmplitude`; the
  leaf prior in stored units, `scale/(k s_f)`, wherever
  [`ConstantGaussianLeaf`](../../src/bartcore/model.hpp) reads `scale/k` (node draws, the branch marginal in
  every move and in grow-from-root's scan, prior draws, `drawLevelShift`, the GIG's leaf precision, and
  the k hyperprior were it reachable); every reader of effective values (`forestTotalFits`, the per-draw
  forest channel, getTrees and printTrees, saved-tree flattening, per-forest predict,
  `forestCalibration`); and the state. Writing `s_f mu` and restoring at `s_f = 1` keeps the wire format;
  a stored `s_f` would be a new optional block under the
  [`stateFormatVersion`](../../src/R_interface_bartcore.cpp) registry rule. `s_f` is itself a
  multiplicative random walk, so it must be folded back into the leaves when it drifts too far, and that
  fold is a bulk transform that needs (a) anyway. Worse, it does not remove the mechanism. Stored-unit
  rounding is relative to the stored leaves, about eps/(a s_f) in effective units, and what enters the
  residual is `a s_f` times the stored gap. The rescale leaves `a s_f` unchanged and an amplitude redraw
  multiplies it by the redraw ratio, which is the same recursion tip's `a g` follows. Floating point is
  scale-free, so (b) is the defect in other coordinates.
- (c) Recompute every forest periodically, every N sweeps. This bounds the amplification without
  removing it, moves every model's draws, costs every user, and would hide the next defect of this class
  rather than expose it.

Recommend (a) with the debug check. (a) is the smallest diff, removes the mechanism, and enforces a
testable invariant. (b) is not a fix, and (c) only bounds the defect. The useful half of (c) is the debug
check. No measured cost changes this; cost only chooses among exact gathers (D3).

D2, the scope. Rule: a cache's gap from its derivation may cross a sweep boundary multiplied never.
(i) A multiplicative transform of leaf values re-derives the cache before the sweep ends. (ii) An
additive transform may update the cache in place: its increment is additive and rounding-level, the same
class as the sweep's own difference updates. (iii) A cache that is rebuilt from scratch every sweep may be
scaled in place between rebuilds; the variance surface and the recorded test fits are such caches.
Recommend applying the fix at the ridge alone. Alternative: also re-derive after the multinomial level
move and the level fibre. That costs a gather over all K forests every multinomial sweep, re-records
multinomial-equivalence, moves every level-Gibbs fit, and buys nothing measured. Evidence that would
change this: the debug check, or the tightened fuzz bound, tripping on a multinomial or level-fibre run.

D3, cost acceptance. Correctness is not traded: the fix lands whatever it costs, and there is no cheaper
correct alternative to fall back to. Measure the plain gather and blocked gathers at widths 4, 8 and 16
(and a row-major pass over tree blocks if the blocked form stalls), and ship the fastest; all give
identical values. Expect 5 to 10 percent of median ms per sweep in the Step 6 BCF cells. Whatever the
number, record it in the landing note and re-check NEWS's "17% for BCF" fused-pass figure. Above 10
percent, report it to VD before landing, as a cost to know about, not a decision to reverse the fix.

D4, what the evidence may reopen. Recommend Steps 7 and 8 as landing evidence, with the pre-registered
criteria applied unchanged:
- The latent SBC arms are admitted only if their ladders now meet
  [Decision 1 - the SBC arms](bcf-latent-evidence.md#decision-1---the-sbc-arms)'s admission clause. The
  sbc.yaml matrix edit is then a follow-up commit.
- The exact gate's statistic (longer batches, AR(1) inflation, seed-spread floor) stays as landed; it can
  only widen. Tightening it again is a follow-up, decided on Step 7's seed-spread measurement.

Two alternatives are rejected. Lifting the exclusion because the gate passes is wrong: the admission
clause is about chain length at large `|a|`, which the gate does not measure. Keeping the exclusion
without re-running is wrong too: the recorded finding may be this defect, and the doc would carry a false
cause.

D5, NEWS. The amplitude family is new in 1.0-0 (main, 0.9-34, has no amplitude code), so no release
carried the defect. But the 1.0-0 section already records fixes only development builds saw: the BCF
zero-multiplier bullet, and the fused-pass entry's note that seeded fits no longer reproduce "an earlier
1.0-0 development build". Recommend a bullet in the same register: amplitude-coupled fits (bcf, and
formula-route multi-forest fits) no longer carry a compounding bias from the ridge move, and their seeded
draws change. Alternative: no bullet, only the D3 re-check; that is inconsistent with the section's own
precedent, and bartCause's 1.0 branch users did run the biased sampler.

## Constraints

- No change to dbarts.h, the state format or rng consumption: the GIG draw sequence is unchanged, and
  draws move through values alone.
- Bitwise gates hold for every configuration with no moved forest (see rng above). The multinomial and
  level-fibre paths are untouched under D2's recommendation.
- Out of scope: the b-move (ridgeB stays off), the gate-statistic and SBC-matrix edits (D4 follow-ups),
  and the variance re-anchor (rule iii).

## Steps

1. Pins, written before the engine change, with their failure on tip recorded (tests/cpp):
   - (a) After a burn-in, write one forest-0 leaf through
     [`Chain::muByTreeForTesting`](../../src/bartcore/chain.hpp) to inject a known gap, then fire the
     ridge. Require `totalFits == gather` bitwise. On tip the gap comes out multiplied by c.
   - (b) Run probit, logistic and gaussian amplitude samplers at the exact gate's shape (one tree per
     forest, one cut) and at an ensemble shape (n = 200, 50 and 25 trees), over several seeds. After
     every sweep, require moved forests to match the gather bitwise. Forests not moved carry additive
     rounding that scales with the forest's working response (residual over its multiplier), not with
     `|total|`, and grows about as the square root of the sweep count: with the fix in, the ridge-off
     treatment forest reached 1.6e4 to 2.6e5 times `m eps (1 + |total|)` at 20000 sweeps. Bound them by
     `C eps max_i |forestY_i| sqrt(sweeps)`, with C set from the fix's measured worst case times 10, and
     record that measurement in the test. On tip, record each seed's worst moved-forest gap and the first
     sweep above 1e-9. Trim the sweep count to whatever still fails on tip by orders of magnitude, within
     about 10 s in total.
   - (c) A `getState`/`setState` round trip at a sweep boundary reproduces forest 0's `totalFits`
     bitwise.
2. Engine, one commit, installed with `--preclean`: the `Forest` gather helper; `rescaleAmplitudeRidge`
   drops the `totalFits` multiply; `afterCombine` gathers every forest it moves; the debug check in
   `Chain::run`. Comments: `Forest` states the rule; the rescale's list of quantities it must keep
   consistent drops `totalFits`; [`ForestCombiner::afterCombine`](../../src/bartcore/combiner.hpp) says
   that an override writing leaves owns rule (i).
3. Update the existing tests/cpp tests:
   - [`testBCFInterweave`](../../tests/cpp/test_sampler.cpp): its 1e-9 sum check becomes bitwise.
   - [`testGeneralAmplitudeRidge`](../../tests/cpp/test_sampler.cpp) builds forests by hand with no
     `leafOf` and a `totalFits` the leaves do not sum to. Give each forest an all-root map and a
     `totalFits` derived from its leaves, then keep the arm that checks the combined fit is invariant.
   - [`fuzzInvariantViolation`](../../tests/cpp/test_fuzz.cpp): check every row rather than every 17th;
     moved forests bitwise, others at Step 1's bound (a `64 m eps (1 + |total|)` bound passes today only at
     29.9 of 64, so it is not used).
   - [`testBCFCombinerSeam`](../../tests/cpp/test_sampler.cpp): its arms build forests by hand with no
     `leafOf`, and its skip cases ("a single occupied leaf", "an all-zero leaf sum") are `update && ridge`
     and assert that `totalFits` is unchanged. Gathering on those no-op returns segfaults there. Give the
     forests maps and a derived `totalFits`, and restate the contract: a no-op return rewrites `totalFits`
     to the gather.
   - Re-run [`testBCFInterweaveKeepTrees`](../../tests/cpp/test_sampler.cpp) and
     [`testBCFGrowForestFromRoot`](../../tests/cpp/test_sampler.cpp).
   - Grep tests/cpp for any other hand-built forest that reaches `afterCombine`.
4. tinytest: the state round trip in ["bartcoreForestFits(restored, 0L)"](../../inst/tinytest/test-bcf.R)
   compares at `tolerance = 1e-5`. Forest 0 becomes `expect_identical` and forest 1 uses 1e-12. Then run
   the whole suite. No file should need a snapshot replay: the four seeded-drift files are single-forest,
   and a grep of the BCF, multinomial and K-forest files finds no literal that depends on draws.
5. Prove the tests discriminate: put the old multiply back in place of the gather. Step 1's pins and the
   debug check must fail. Revert and `touch` the file.
6. Cost, run by the maintainer on a quiet machine, arm64 and the x86 box, never under other load:
   - Add BCF cells to [`runScenarios`](../../benchmarks/R/bench-sampler.R): gaussian and probit bcf, n =
     1e3 and 1e4, p = 10, 200 prognostic and 50 treatment trees, one chain, one thread.
   - Record on tip and on the slice from two private libraries, alternating ABAB for five pairs. Report
     the median ms per sweep.
   - `bench-sampler.R compare` on the existing cells must sit within noise; their code does not change.
7. Exact gates: every gate in exact-gates.yaml's list in quick mode.
   [`exactLatentBCF`](../../benchmarks/R/bcf-latent-exact.R) in quick and in full mode (local, about 40
   minutes), and bcf-exact.R, bcf-exact-weak.R and bcf-exact-restricted.R in full mode.
   - Oracle for the re-records: a `pooled` mode for bcf-latent-exact.R that runs the quick shape over at
     least 300 seeds and scores each gated channel's pooled mean against the quadrature, using the seed
     spread as its se. Require `|z| < 3` on every channel (7.1 on E[tau] at tip).
   - Measure mode 2a's 20-seed spread against its batch se again (recorded at 8x to 30x understatement).
8. Calibration: the burn-bcf-probit and burn-bcf-logistic ladders; the bcf-probit and bcf-logistic
   points at R = 200 with their n = 40 controls; and the gaussian bcf arm at its recorded settings
   ([`sbcBurnSweeps`](../../benchmarks/R/sbc.R)). Record the verdicts; D4 governs admission.
9. Equivalence. The stored baselines are host-bound: on arm64, tip against `equivalence-d2b9827a` already
   moves five scenarios (friedman, probit, weighted, splitprobs, quants, |z| up to 2.59). So compare the
   fix against baselines recorded on tip on the same host, then re-record the stored ones on their
   recording host.
   - equivalence.R `--strict-coverage`: every scenario identical except bart2twoforest.
     multinomial-equivalence: 11 identical. bcf-equivalence: every scenario that draws `a` moves, which
     is all 15 at their defaults; confirm each.
   - Every `|z|` flag must fall on amplitude-coupled channels, and Step 7 adjudicates it.
   - Re-record equivalence.R and bcf-equivalence. Their MANIFEST rows name Step 7 as the oracle (P17).
     The same push updates the baseline filenames in cpp-tests.yaml, equivalence.yaml and exact-gates.yaml
     (the bcf cross-host compare included).
10. Sanitizers: tests/cpp under ASAN and UBSAN, per [Gate hygiene](README.md#gate-hygiene).
11. Records, below.

## Verification

- `cd tests/cpp && make && ./test_bartcore`: all pass, including Step 1. After Step 5's mutation, Step 1
  fails.
- `R CMD INSTALL --preclean -l <lib> .`, then `R_LIBS=<lib> Rscript -e
  'tinytest::test_package("dbarts")'`: no failures.
- `Rscript benchmarks/R/bcf-latent-exact.R quick` and `Rscript benchmarks/R/bcf-latent-exact.R`: every
  matched quantity within the gate's `|z| <= 4`. The pooled mode shows `|z| < 3` on every channel.
- The debug check: tests/cpp under `!NDEBUG` reports no trip.
- `Rscript benchmarks/R/equivalence.R compare <tip baseline, same host> --strict-coverage`: every
  scenario "identical draws (same RNG stream)" except bart2twoforest.
  `multinomial-equivalence.R compare <tip baseline, same host>`: 11 identical. After the re-records, every compare
  reports its full scenario count identical.
- `Rscript benchmarks/R/bench-sampler.R compare benchmarks/baselines/bench-sampler-127f04ee.csv`: exit 0.
  The Step 6 BCF A/B within D3.

## Records at landing

- [The ASIS ridge](../design/multiplier-combiner.md#the-asis-ridge): the design note for this
  posterior-changing change. It states the rule, drops `totalFits` from the list of quantities the move
  keeps consistent, and bumps the Status line.
- docs/design/bcf.md: [Burn-in under strong prognostic signal (2026-07-10)](../design/bcf.md#burn-in-under-strong-prognostic-signal-2026-07-10)
  and [Calibration (2026-07-07)](../design/bcf.md#calibration-2026-07-07) are revised wherever Step 8
  changes their reading. [Exact-posterior gate](../design/bcf.md#exact-posterior-gate) gets the Step 7
  result.
- docs/plans/bcf-latent-evidence.md: its "metastable" explanation of mode 2a (under Decision 2) becomes
  what Step 7 measured, with the drift named as its cause if the spread closes. Decision 1's finding is
  amended with Step 8. The comment in ["metastable"](../../benchmarks/R/bcf-latent-exact.R) and the
  ["Not in the matrix: BCF"](../../.github/workflows/sbc.yaml) exclusion note change to match.
- TODO: forest-cache-drift is removed. bcf-sigma-tail-mixing is re-read against the gaussian arm's
  Step 8 result.
- NEWS per D5. This file's Status line and Landing note, and its INDEX row.
