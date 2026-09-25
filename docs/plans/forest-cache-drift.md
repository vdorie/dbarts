# forest-cache-drift

Status: IN PROGRESS, 2026-09-24. The removal (dec-B127 in [docs/decisions.md](../decisions.md)) is implemented
(ac7f70b7, see the Landing note); D4 and D5 below are open for VD, and the gaussian equivalence baseline is owed from
its recording host.
agent: opus (engine, component tests, gate battery, re-records); sonnet (records: docs, INDEX, TODO)
rng: posterior-changing. The amplitude rescaling move and its per-sweep GIG draw are removed, so
  every chain with an updating scale-mixture amplitude (bcf's prognostic forest, and each
  formula-route forest without a basis) draws differently from its first sweep on, and the
  implemented chain's stationary law moves, because the defect biased it. Every other
  configuration (single-forest, heteroscedastic, multinomial, and an amplitude coupling whose
  scale-mixture forests are all held) must replay its baseline bitwise.
window: none (VD 2026-09-24: the right fix, whatever its timing against the merge).
budget: engine mostly deletion (combiner.hpp, chain.hpp, the bridge, random.c) plus a ~30-line debug
  check; tests/cpp ~250 (two pins, tests rewritten or removed); tinytest ~30; benchmarks ~60 (a
  pooled-seed gate mode); two baseline re-records; records.

## Goal

The amplitude rescaling move is gone, and with it the only transform that multiplied a forest's
cached fits separately from its leaves. Every forest cache differs from its derivation by
accumulated additive rounding only, a state restore at a sweep boundary reproduces the cache to that
rounding, and a debug-build check enforces the rule for any later transform. The latent BCF exact
gate passes in both modes, and the gate and calibration records that read the drift as mixing are
measured again.

## Context

The defect. retired: [`AmplitudeForestCombiner::rescaleAmplitudeRidge`](../../src/bartcore/combiner.hpp),
now removed, moved forest f along its ridge, `(a_f, leaves) -> (a_f/c, c leaves)`, once per sweep,
from [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp)'s `afterCombine`. It multiplied
`muByTree` by c, and separately multiplied the [`Forest`](../../src/bartcore/combiner.hpp)'s
`totalFits` by c. Let `g_i = totalFits_i - sum_t mu_t[leafOf_t(i)]`. The sweep keeps `totalFits` by
difference updates ([`Chain::rollTreeResidual`](../../src/bartcore/chain.hpp),
[`Chain::fusedRollPass`](../../src/bartcore/chain.hpp),
[`Chain::finalizeTotalFits`](../../src/bartcore/chain.hpp)), which carry g unchanged through every leaf
redraw and add a rounding increment. The rescale mapped g to `c g`. `a_f` is redrawn from its
conditional between rescales, so the c factors are not ratios of one quantity and do not cancel: g
followed a multiplicative random walk while the fits stayed O(1).

This had two effects. First, g enters every tree's partial residual (tree t fits `y - g - sum_{s != t}
f_s`), so the leaves absorb `-g` while `totalFits` still tracks the data. The sigma draw (through the
combined location), the latent refresh, the amplitude conditional and the per-forest reporting
channel then read one function (the cache). The leaf prior, the move's own GIG statistic, the saved
trees, predict and the state read another (the leaves). Second,
[`Chain::setState`](../../src/bartcore/chain.hpp) and [`Chain::installForest`](../../src/bartcore/chain.hpp)
re-derive the cache ([`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp)). So a restored chain
dropped g without warning, and continued from a different state than the one that was saved.

Measured by the investigation's probes. At the exact gate's single-tree shape, the index-scale gap
`|a| max|fit - leaf|` grew from 1e-11 to 0.3 under probit (seed 13). Within 20000 sweeps it reached
0.53 under logistic and 2.4e-3 under gaussian. Latent BCF tau was biased +0.005 at the gate's quick
length and +0.027 at full length, a pooled z of 7.1 over 303 seeds. That bias is what
[`exactLatentBCF`](../../benchmarks/R/bcf-latent-exact.R) caught in quick mode on arm64, and it does
not depend on the host. A scratch patch that rebuilt `totalFits` from the rescaled leaves zeroed the
gap and restored the gate.

Which models were affected: every forest that the bridge gave a half-Cauchy amplitude,
[`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp) having switched the move on exactly when
`amplitudePriorScale > 0`. That covers bcf's prognostic forest and, in the formula route's K-forest
fits, each forest without a basis (`forestParams` in [R/model.R](../../R/model.R); in bart2twoforest,
forest 0 only), under gaussian, probit and logistic.

Design: [The ASIS ridge](../design/multiplier-combiner.md#the-asis-ridge),
[5. Correctness](../design/level-fibre.md#5-correctness),
[Identification and the level-centering move](../design/multinomial.md#identification-and-the-level-centering-move),
[Decision 2 - the exact gate](bcf-latent-evidence.md#decision-2---the-exact-gate).

## Audit

This covers every site in src/bartcore and the bridge that transforms leaf values, a calibration, a
response scale or a forest's level in bulk, and what each does to derived caches. "Amplifies" means
the gap is multiplied. "Preserves" means the gap is carried with additive rounding. "Clears" means
the cache is re-derived from scratch.

| site | transform | caches, as adjusted | gap | action |
|---|---|---|---|---|
| the rescaling move (removed) | leaves times c | `totalFits` times c, separately; the recorded test fits and the keepTrees slot times c | amplified, never cleared | removed |
| [`MultinomialForestCombiner::afterCombine`](../../src/bartcore/combiner.hpp) | occupied leaves plus c/m | `totalFits` plus c (m fl(c/m) is not exactly c); test fits untouched by design | preserves: additive, rounding-level | none (rule ii) |
| [`Chain::drawLevelShift`](../../src/bartcore/chain.hpp) | occupied leaves plus c_t, where the c_t sum to 0 | `totalFits` untouched, correct algebraically | preserves: the projection's rounding residual | none (rule ii) |
| [`Chain::reanchorVarianceForest`](../../src/bartcore/chain.hpp) | variance factors times g | `combinedVariance` and its test twin times f, where g^m' is not exactly f | preserves for one sweep; [`Chain::sweepVarianceForest`](../../src/bartcore/chain.hpp)'s closing product and [`Chain::refreshVarianceTestFits`](../../src/bartcore/chain.hpp) clear it | none (rule iii) |
| [`GaussianResponse::setResponse`](../../src/bartcore/model.hpp), [`GaussianResponse::setOffset`](../../src/bartcore/model.hpp), [`GaussianResponse::setData`](../../src/bartcore/model.hpp), [`GaussianResponse::restoreScale`](../../src/bartcore/model.hpp) | the response transform moves; leaves are reinterpreted, not transformed | none | none | none |
| [`Chain::setForestPriorScale`](../../src/bartcore/chain.hpp), [`Chain::setModel`](../../src/bartcore/chain.hpp), [`Chain::setForestBasis`](../../src/bartcore/chain.hpp) | leaf scale, k | none: no drawn leaf is touched | none | none |
| [`Chain::rebuildFitsFromParameters`](../../src/bartcore/chain.hpp), and the per-observation update session, which ends in the same rebuild | per-tree re-route, leaves unchanged | `totalFits` subtract, then add | preserves | none |
| sweep upkeep: the three roll functions above, [`VarianceForest::applyLeafFactor`](../../src/bartcore/chain.hpp) | per-tree redraws | difference updates | preserves (1.2e-15 at 200 trees, [`runEnsembleTests`](../../tests/cpp/test_ensemble.cpp)) | none |
| [`Chain::sampleTreesFromPrior`](../../src/bartcore/chain.hpp), [`Chain::sampleNodeParametersFromPrior`](../../src/bartcore/chain.hpp), [`Chain::applyNewData`](../../src/bartcore/chain.hpp), [`Chain::forceRefreshTrees`](../../src/bartcore/chain.hpp), [`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp), [`Chain::rebuildLiveForestRemapped`](../../src/bartcore/chain.hpp), [`Chain::sampleVarianceForestFromPrior`](../../src/bartcore/chain.hpp), [`Chain::rebuildVarianceForest`](../../src/bartcore/chain.hpp), [`Chain::refreshVarianceForest`](../../src/bartcore/chain.hpp) | wholesale replacement | rebuilt from zero | clears | none |
| the bridge, [`applyAmplitudeSpec`](../../src/R_interface_bartcore.cpp) and the state entries | none: forest state is only passed through, and state goes through the re-deriving restore | reporting reads the cache (`$getForestFits`) and the leaves (trees, predict) as they are | none of its own | none |

The remaining caches are rebuilt every sweep or scratch, and no transform reaches them: node
sufficient statistics and `treeY` (rebuilt by the roll, dead across sweep boundaries), the
multinomial offset slab (rematerialized at sweep entry and at reporting), the linear leaf's U'WU
cache (no bulk transform of vector leaves exists), and combiner scratch.

The rule, which stands with the move gone: a cache's gap from its derivation may never cross a
sweep boundary multiplied. (i) A multiplicative transform of leaf values re-derives the cache
before the sweep ends. (ii) An additive transform may update the cache in place: its increment is
additive and rounding-level, the same class as the sweep's own difference updates. (iii) A cache
that is rebuilt from scratch every sweep may be scaled in place between rebuilds; the variance
surface and the recorded test fits are such caches. No multiplicative leaf transform remains, so
rule (i) binds only a future one, a restored rescaling move included. The TODO item
level-shift-capi would add the first transform a host drives; it is additive, so rule (ii) covers
it.

## Decision

The move is removed (dec-B127). With the defect fixed by re-deriving the cache after every
rescale, a comparison against the same model without the move, over seven designs at two sample
sizes and 20 seeds each, gave at most 1.30 times the effective draws per second on the residual
scale or the treatment effect (gaussian, strong signal, n = 1000, treatment effect; bootstrap
interval 0.83 to 1.70), against a bar of 1.5 set before the run. The move with its fix cost 4 to 6
percent of every sweep. The amplitude itself mixed better in some designs (up to 1.56 times, in
the formula-route design), but no quantity a user reads did; without the move the amplitude still
reaches its stationary level within a few hundred sweeps in the strong-signal large-amplitude
design, against a residual-scale transient of tens of thousands on both builds. The slow mixing in
the strong-signal designs is the tree-structure limit, which the move does not touch. The two
builds' posterior means agreed (largest paired t over 88 design-channel cells 2.23, about what
chance gives).

The move stays on record as an option for a later mixing experiment:
[The ASIS ridge](../design/multiplier-combiner.md#the-asis-ridge) keeps its derivation, and the
history keeps its code. A restoration must re-derive the moved forest's cached fits from its leaves
after every rescale and never multiply a cache separately.

Two questions are open for VD.

D4, what the evidence may reopen. Recommend Steps 6 and 7 as landing evidence, with the
pre-registered criteria applied unchanged:
- The latent SBC arms are admitted only if their ladders now meet
  [Decision 1 - the SBC arms](bcf-latent-evidence.md#decision-1---the-sbc-arms)'s admission clause.
  The sbc.yaml matrix edit is then a follow-up commit.
- The exact gate's statistic (longer batches, AR(1) inflation, seed-spread floor) stays as landed;
  it can only widen. Tightening it again is a follow-up, decided on Step 6's seed-spread
  measurement.

Two alternatives are rejected. Lifting the exclusion because the gate passes is wrong: the
admission clause is about chain length at large `|a|`, which the gate does not measure. Keeping the
exclusion without re-running is wrong too: the recorded finding may be this defect, and the doc
would carry a false cause.

D5, NEWS. The amplitude family is new in 1.0-0 (main, 0.9-34, has no amplitude code), so no release
carried the defect or the move. But the 1.0-0 section already records changes only development
builds saw: the BCF zero-multiplier bullet, and the fused-pass entry's note that seeded fits no
longer reproduce "an earlier 1.0-0 development build". Recommend a bullet in the same register:
amplitude-coupled fits (bcf, and formula-route multi-forest fits) no longer carry a compounding bias
from the amplitude rescaling move, which is removed; their seeded draws change and each sweep is
about 5 percent faster. Alternative: no bullet; that is inconsistent with the section's own
precedent, and bartCause's 1.0 branch users did run the biased sampler.

## Constraints

- No change to dbarts.h or the state format. The amplitude and its half-Cauchy auxiliary are still
  drawn from their exact conditionals every sweep.
- Bitwise gates hold for every configuration with no updating scale-mixture amplitude (see rng
  above). The multinomial and level-fibre paths are untouched.
- Out of scope: the gate-statistic and SBC-matrix edits (D4 follow-ups), and the variance
  re-anchor (rule iii).

## Steps

1. Pins, written before the engine change, with their failure on tip recorded (tests/cpp,
   [`testAmplitudeCacheDrift`](../../tests/cpp/test_sampler.cpp),
   [`testAmplitudeCacheRestore`](../../tests/cpp/test_sampler.cpp)):
   - (a) Drift: probit, logistic and gaussian amplitude samplers at the exact gate's shape (one
     tree per forest) and at an ensemble shape (n = 200, 50 and 25 trees), over several seeds.
     After every sweep, every forest's `totalFits` is within `C eps max_i |forestY_i| sqrt(sweeps)`
     of its leaf gather, with C set from the measured worst case times 10 and that measurement
     recorded in the test. On tip, record each configuration's worst gap and the first sweep above
     the bound.
   - (b) Restore: a `getState`/`setState` round trip at a sweep boundary on an amplitude-coupled
     sampler reproduces each forest's `totalFits` to an additive-rounding bound.
2. Engine, one commit, installed with `--preclean`: remove the move, the `ridge` fields and every
   branch that reads them, the keepTrees slot and test-fit rescale, the testing hook that fired
   the move, the bridge's derivation of the flag, and the GIG generator nothing else uses. Add the
   debug check to [`Chain::run`](../../src/bartcore/chain.hpp): at the end of every sweep under
   `!NDEBUG`, every constant-leaf forest of a double-precision chain within `1e-8 (1 +
   |totalFits_i|)` of its gather. Float chains are exempt: the fp32 tests carry gaps above 1e-10 by
   design. The R build defines NDEBUG and compiles the check out.
3. Tests: remove or rewrite every test that exists for the move (the GIG moment case, the
   interweave tests, the ridge arms of the combiner seam and grow-from-root pins, the general
   ridge test, the tinytest ridge blocks), keep coverage of the amplitude draw, check every row in
   [`fuzzInvariantViolation`](../../tests/cpp/test_fuzz.cpp) at the Step 1 bound, and tighten the
   state round trip in ["bartcoreForestFits(restored, 0L)"](../../inst/tinytest/test-bcf.R) to
   1e-12.
4. Prove the pins discriminate: on tip, Step 1 (a) fails.
5. Sanitizers: tests/cpp under ASAN and UBSAN, per [Gate hygiene](README.md#gate-hygiene).
6. Exact gates: every gate in exact-gates.yaml's list in quick mode.
   [`exactLatentBCF`](../../benchmarks/R/bcf-latent-exact.R) in quick and in full mode, and
   bcf-exact.R, bcf-exact-weak.R and bcf-exact-restricted.R in full mode.
   - Oracle for the re-records: a `pooled` mode for bcf-latent-exact.R that runs the quick shape
     over at least 300 seeds and scores each gated channel's pooled mean against the quadrature,
     using the seed spread as its se. Require `|z| < 3` on every channel (7.1 on E[tau] at tip).
   - Measure mode 2a's 20-seed spread against its batch se again (recorded at 8x to 30x
     understatement).
7. Calibration: the burn-bcf-probit and burn-bcf-logistic ladders; the bcf-probit and bcf-logistic
   points at R = 200 with their n = 40 controls; and the gaussian bcf arm at its recorded settings
   ([`sbcBurnSweeps`](../../benchmarks/R/sbc.R)). Record the verdicts; D4 governs admission.
8. Equivalence. The stored baselines are host-bound: on arm64, tip against `equivalence-d2b9827a`
   already moves five scenarios. So compare the slice against baselines recorded on tip on the same
   host, then re-record the stored ones under the MANIFEST's convention.
   - equivalence.R `--strict-coverage`: every scenario identical except bart2twoforest.
     multinomial-equivalence: 11 identical. bcf-equivalence: every scenario that draws `a` moves,
     which is all 15 at their defaults; confirm each.
   - Every `|z|` flag must fall on amplitude-coupled channels, and Step 6 adjudicates it.
   - Re-record equivalence.R and bcf-equivalence. Their MANIFEST rows name Step 6 as the oracle
     (P17). The same push updates the baseline filenames in cpp-tests.yaml, equivalence.yaml and
     exact-gates.yaml (the bcf cross-host compare included).
9. Speed: `bench-sampler.R compare` against the current baseline; its cells are single-forest and
   do not change.
10. Records, below.

## Verification

- `cd tests/cpp && make && ./test_bartcore`: all pass, including Step 1, with the debug check live
  (tests/cpp does not define NDEBUG). On tip, Step 1 (a) fails.
- `R CMD INSTALL --preclean -l <lib> .`, then `R_LIBS=<lib> Rscript -e
  'tinytest::test_package("dbarts")'`: no failures.
- `Rscript benchmarks/R/bcf-latent-exact.R quick` and `Rscript benchmarks/R/bcf-latent-exact.R`:
  every matched quantity within the gate's `|z| <= 4`. The pooled mode shows `|z| < 3` on every
  channel.
- `Rscript benchmarks/R/equivalence.R compare <tip baseline, same host> --strict-coverage`: every
  scenario "identical draws (same RNG stream)" except bart2twoforest.
  `multinomial-equivalence.R compare <tip baseline, same host>`: 11 identical. After the re-records,
  every compare reports its full scenario count identical.
- `Rscript benchmarks/R/bench-sampler.R compare benchmarks/baselines/bench-sampler-127f04ee.csv`:
  exit 0.

## Records at landing

- [The ASIS ridge](../design/multiplier-combiner.md#the-asis-ridge) becomes the record of a removed
  move: what it was, why it was removed, the last commit that carried it, the rule a restoration
  must follow, and the exponent-rule derivation.
- docs/design/bcf.md: wherever it describes the move or its mixing effect, and a pointer in its
  mixing discussion that the move is an available option. [Exact-posterior gate](../design/bcf.md#exact-posterior-gate)
  gets the Step 6 result.
- docs/plans/bcf-latent-evidence.md: its "metastable" explanation of mode 2a (under Decision 2)
  becomes what Step 6 measured. Decision 1's finding is amended with Step 7. The comment on
  [`batchStats`](../../benchmarks/R/bcf-latent-exact.R) and the
  ["Not in the matrix: BCF"](../../.github/workflows/sbc.yaml) exclusion note change to match.
- docs/plans/bartcore-landing: rows that describe the move.
- docs/architecture.md, benchmarks/README.md, the feature matrix if it lists the move.
- TODO: forest-cache-drift is removed. bcf-sigma-tail-mixing is re-read against the gaussian arm's
  Step 7 result.
- NEWS per D5. This file's Status line and Landing note, and its INDEX row.


## Landing note (2026-09-24)

The engine (ac7f70b7), the harness (af402611), the BCF baseline (90f9bccc) and these records; D4 and D5 stay open, so
the Status stays IN PROGRESS.

Engine, ac7f70b7. The move is gone from [`AmplitudeForestCombiner`](../../src/bartcore/combiner.hpp), which no
longer overrides `afterCombine`; with it went the per-forest `ridge` flags and the bridge's derivation of them, the
keepTrees slot and test-fit rescale, the testing hook that fired the move, and the GIG generator. `afterCombine` returns
nothing now, its value having been read only by the move's tests, and its Doxygen states the forest cache rule.
[`Chain::run`](../../src/bartcore/chain.hpp) checks the rule under `!NDEBUG`, which tests/cpp compiles with; restoring
the move verbatim trips it. Every amplitude-coupled chain draws differently from its first sweep; dbarts.h and the state
format are unchanged.

Pins. [`testAmplitudeCacheDrift`](../../tests/cpp/test_sampler.cpp) and
[`testAmplitudeCacheRestore`](../../tests/cpp/test_sampler.cpp) measure each forest's gap in units of
`eps max_s max_i |forestY_i(s)| sqrt(sweeps)`, the running maximum of the working response over the sweeps seen, since
a gap keeps the rounding of a sweep whose multiplier was small. C = 2000: ten times the worst measured without the move,
184, from 1000 fuzz seeds (the pins' own worst is 5.0). On the previous tip the drift pin fails at the latent gate's
shape, probit at ratio 4.7e10 (first over the bound at sweep 1924) and logistic at 7.7e11 (sweep 2084), and the restore
pin at 1.12e4 (probit, gate shape); gaussian and the ensemble shape stay under the bound there within the pins' sweep
counts. [`fuzzInvariantViolation`](../../tests/cpp/test_fuzz.cpp) checks every row at the same bound, its burn-in run a
sweep at a time. tests/cpp passes, and under ASAN and UBSAN with no diagnostic.

Harness, af402611. The latent BCF SBC arms build their host under their own link; before it none of them ran. The
latent gate's `pooled` run reports each channel's seed spread over its batch se.

Gates. Every exact-gates.yaml gate passes in quick mode. bcf-exact.R, bcf-exact-weak.R and bcf-exact-restricted.R pass
in full mode, bcf-latent-exact.R in quick and full mode, and its pooled run at 300 seeds sits at worst `|z|` 0.63
(probit) and 1.12 (logistic), 7.1 at the previous tip. Mode 2a's seed spread is 0.7x to 1.2x its batch se over 20
seeds; [Decision 2 - the exact gate](bcf-latent-evidence.md#decision-2---the-exact-gate) records what that retires.
The full tinytest suite passes, 8887 tests, with no snapshot replayed.

Equivalence, against baselines recorded on the previous tip on the same host: equivalence.R 52 of 53 identical under
`--strict-coverage`, bart2twoforest moving at max `|z|` 1.00; multinomial 11 of 11 identical; bcf all 15 moving, every
flag on an amplitude-coupled channel. bcf-equivalence is re-recorded as `bcf-equivalence-ac7f70b7.rds` (90f9bccc), the
exact gates its oracle; this host reproduces the stored BCF baseline bitwise at the previous tip. equivalence.R is NOT re-recorded:
this host reproduces `equivalence-d2b9827a.rds` in 48 of 53 scenarios only, so a recording here would move five
scenarios CI compares bitwise. It is owed from the stored baselines' recording host, and until then cpp-tests.yaml's
bitwise compare fails on bart2twoforest.

Calibration. See [Decision 1 - the SBC arms](bcf-latent-evidence.md#decision-1---the-sbc-arms)'s re-measurement: the
ladders still fail the admission clause at `|a| >= 5` on both links, with and without the move; the `R = 200` verdicts
are 10 of 13 (probit) and 12 of 13 (logistic), every functional inside the matrix band, and the `n = 40` controls pass
13 of 13. The gaussian arm at its recorded settings passes 13 of 15, sigma included, against 9 of 15 with the move on
the same host ([Calibration (2026-07-07)](../design/bcf.md#calibration-2026-07-07)).

Speed. Not measured: the host was loaded (1-minute load 6 to 7 from other work). `bench-sampler.R compare` against
`bench-sampler-127f04ee.csv`, alternating three rounds on each build, flagged zero to four cells at 1.05 to 1.11 on
the previous tip and one to four on the slice, `embedded-offset-run1-n1000-t75` the most often on both. The single-forest
cells' code does not change, the check being compiled out under NDEBUG; the compare is owed on a quiet machine.
