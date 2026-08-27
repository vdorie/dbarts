# gate-hardening-1.0

agent: one implementer (worktree wt/gate-hardening-1.0)
rng: neutral (gates only; no engine or package code changes; existing 18
     equivalence scenarios stay bit-identical against equivalence-cf99a00)
budget: new gates for the poison-sweep blind spots (gate-blindspot-audit.md)

Deliverable of review 3 (gate-blindspot-audit.md): each sub-item exists
because a specific deliberate kernel breakage slipped every R-side gate
or all gates.

## Defect list (poison numbers from the audit matrix)

1. Poison 13 - BCF a-glue prior precision (chain.hpp:2049 drops
   1/aVariance): TRUE survivor, zero gates. DESIGN ONLY this turn (an
   SBC experiment, bcf-weak, decides between the two designs below).
2. Poison 16 - GP nugget drops /w_i (model.hpp:721): caught only by
   tests/cpp; the gp equivalence scenario is unweighted. Gate: new
   `wtgp` equivalence scenario (GP leaves x non-unit weights).
3. Poison 12 - grouped precision count (model.hpp, drawGroupEffects
   weight accumulation `weightScratch[j] += w` -> `+= 1.0`): cpp-only;
   tinytest's "catch" was a CI hang. Gate: new `grouped` equivalence
   scenario (rbart_vi in-core, gamma tau prior, non-unit weights).
4. Poison 15 - linear ridge drops sigma^2 (model.hpp:304, score side
   only): no exact gate. Gate: benchmarks/R/linear-exact.R.
5. Poison 2 - birth/death reverse node-selection count on the wrong
   tree (moves.hpp ~245): slipped cpp AND change-balance. Gate:
   benchmarks/R/bd-balance.R (exact-enumeration birth/death
   detailed-balance arm).
6. Poison 8 - chi-k shape mislabel (model.hpp:1729, shape
   0.5(M + 2nu - 1) instead of 0.5(M + nu)): cpp-only; chik's nu = 1.5
   shifts the shape by only +0.5. Gate: new `chik2` equivalence
   scenario with a df chosen so the shape error is separable.

## Sub-item 1 design (poison 13, BCF a-glue prior precision)

UNBLOCKED 2026-07-09: the bcf-weak SBC found the prior-dominant BCF
configuration (n = 40) well-behaved (sigma dead-centre; |a| channel
ecdf 0.043), so the prior-dominated design is implemented as
benchmarks/R/bcf-exact-weak.R (standalone companion to bcf-exact.R,
which stays untouched):

- Small n (~40), K = 2 predictor cells, weak mu signal, a free
  (Cauchy(0, sd.control) via the IG scale mixture), b fixed (0, 1),
  sd.control chosen SMALL (0.5) so the prior precision 1/aVariance
  dominates the data precision sum_i w mu_i^2 / sigma^2 - the regime
  where poison 13 (dropping the prior precision from the a full
  conditional) inflates the a posterior spread massively instead of
  shifting E[a mu] by 1e-4.
- Observables: E[a mu] per cell (identified mean check) plus
  P(|a| <= q) at fixed thresholds - a SPREAD check, because the
  ridge-regime a posterior has Cauchy tails (E[|a|] is not finite), so
  probabilities/quantiles are the well-posed variance-channel gate the
  mean-only check was blind to.
- Exact side: tree spaces at K = 2 are tiny (stump or one split per
  forest); conditional on (a, sigma) the leaf block is closed form;
  quadrature over sigma and over a through the substitution
  a = sd.control * tan(t) (the Cauchy prior becomes uniform in t), so
  the a integral covers the whole real line with no truncation.
- Chain settings conservative per the SBC caution (the (a, mu) scale
  ridge mixes slowly at larger n): generous burn-in, |a|
  autocorrelation checked at this n before fixing draw counts.

## Working notes (running)

- Private library: .Rlib inside the worktree; ~/.Renviron overrides
  R_LIBS_USER, so all installs/runs use R_LIBS=<worktree>/.Rlib and
  R CMD INSTALL -l <worktree>/.Rlib. Verified dbarts loads from it.
- Clean preclean install: ~25s on this machine.
- equivalence.R compare already tolerates baseline-missing scenarios
  (the `uncovered` path: reports coverage, warns, exit 0 unless
  --strict-coverage). Requirement met by existing code.
- linear-exact derivation validated against the clean sampler
  (max gap 0.0012 at 20k draws, K=4 well-separated cells).
- Poison-15 sensitivity tuning (marginal-side ridge poison only; the
  draw at model.hpp:341 stays clean): the poison mostly reweights tree
  sizes, so the gate checks the leaf-count distribution (exact
  enumeration gives it) as well as E[fit]. Config nPer=10, base=0.8,
  sd=0.6 gives clean size dist {1: .47, 2: .45, 3: .08} and
  E[fit] poison gap ~0.011.

## Landing record

### Sub-item 4: benchmarks/R/linear-exact.R (poison 15)

Gate: single tree, one ordinal predictor, K = 4 cells (n = 40), linear
leaf with fixed sigma; 15 enumerable trees, closed-form ridge marginal
per leaf. Checks P(#leaves) and E[fit] per cell with batch-means z,
|z| < 4. Full mode 200k kept draws, ~2.5 s.

Clean run (full, 200k draws): all |z| <= 1.5, exit 0.
  P(1 leaves) engine 0.4688 exact 0.4674 z +1.2
  P(2 leaves) engine 0.4498 exact 0.4510 z -1.1
  P(3+ leaves) engine 0.0814 exact 0.0816 z -0.2
  E[f(x=1..4)] z +0.9 -1.5 -1.3 -0.3
  OK: linear-leaf sampler matches the exact posterior

Fail-on-poisoned proof (poison 15: model.hpp:304
`ridge = (k/scale)^2 * residualVariance` -> `(k/scale)^2`, score side
only; R CMD INSTALL --preclean into the private .Rlib):
  P(1 leaves)      0.3579     0.4674    0.00225    -48.5 <- FAIL
  P(2 leaves)      0.5086     0.4510    0.00238    +24.2 <- FAIL
  P(3+ leaves)     0.1335     0.0816    0.00160    +32.3 <- FAIL
  E[f(x=1)]       -1.0254    -1.0313    0.00128     +4.6 <- FAIL
  E[f(x=2)]       -0.5118    -0.5017    0.00105     -9.5 <- FAIL
  E[f(x=3)]        0.1208     0.1129    0.00108     +7.3 <- FAIL
  E[f(x=4)]        0.6815     0.6844    0.00127     -2.3
  FAIL: linear-leaf sampler deviates from the exact posterior (exit 1)
Poison reverted (git checkout src/bartcore/model.hpp), --preclean
reinstall, gate re-passes (quick, all |z| <= 0.8, exit 0).

### Sub-item 5: benchmarks/R/bd-balance.R (poison 2)

Gate: single tree, one ordinal predictor, K = 4 cells (n = 40),
constant leaves, fixed sigma, birth/death-dominated kernel
(birth_death = 0.99, change = 0.01, swap = 0; the bridge requires
birth_death < 1 and change targets the same posterior). Fully
enumerable tree space; the exact posterior aggregates onto the 8
cut-set partitions, which the engine arm observes per kept sample from
getTrees split values. |z| < 4 per partition, batch-means MC error.
Full mode 300k kept draws, ~3 s.

Clean run (full, 300k draws): all |z| <= 1.3, exit 0.
  3 0.4179/0.4166 +1.3; (none) 0.1872/0.1875 -0.3; 1+3 0.1635/0.1636
  -0.3; 1 0.1048/0.1054 -0.9; 2 0.0429/0.0433 -1.1; 2+3 0.0394/0.0394
  +0.0; 1+2 0.0258/0.0260 -0.8; 1+2+3 0.0184/0.0181 +1.2
  OK: birth/death chain matches the exact posterior over partitions

Fail-on-poisoned proof (poison 2: moves.hpp death branch, reverse
node-selection count probabilityOfSelectingNodeForBirth computed on
the PRE-death tree, i.e. moved above tree.orphanChildren; --preclean
into the private .Rlib):
  partition      engine      exact       MCse        z
  3              0.0000     0.4166    0.00000     -Inf <- FAIL
  (none)         0.0000     0.1875    0.00000     -Inf <- FAIL
  1+3            0.0000     0.1636    0.00000     -Inf <- FAIL
  1              0.0000     0.1054    0.00000     -Inf <- FAIL
  2              0.0000     0.0433    0.00000     -Inf <- FAIL
  2+3            0.0000     0.0394    0.00000     -Inf <- FAIL
  1+2            0.0000     0.0260    0.00000     -Inf <- FAIL
  1+2+3          1.0000     0.0181    0.00000     +Inf <- FAIL
  FAIL: birth/death stationary distribution deviates from exact (exit 1)
In this single-variable configuration the poison is absorbing: at the
full tree every pre-death leaf is single-cell (not birthable), the
wrong-tree count returns 0, and deaths are never accepted - the gate
z-fails rather than hanging. Poison reverted (git checkout
src/bartcore/moves.hpp), --preclean reinstall, gate re-passes
(quick, all |z| <= 2.1, exit 0).

### Sub-items 2, 3, 6: equivalence scenarios wtgp, grouped, chik2

equivalence.R additions (RNG-neutral; verified by recording a quick
baseline with the HEAD script and comparing with the modified script:
all 18 existing scenarios "identical draws (same RNG stream)", the 3
new ones reported as uncovered with exit 0):
- EQUIVALENCE_SCENARIOS env filter for targeted record/compare passes.
- wtgp: GP leaves x strictly positive non-unit weights (poison 16's
  blind spot; the weighted score path's nugget, not the zero-weight
  fallback).
- grouped: rbart_vi in-core grouped intercepts, gamma tau prior (NOT
  half-Cauchy: SBC found its tail can stall the tau slice sampler),
  weights on [0.5, 3]; new fitViaRbart path; summaries add
  tau.mean/tau.sd/ranef.1-8.
- chik2: normal(chi(50)) at 20 trees (scenario nTrees override) so the
  chi-k shape term is separable.
- compare gains a DISJOINT-SEED-RANGE failure channel alongside the
  Welch z: a diverged run inflates the Welch variance so |z| caps near
  sqrt(n.seeds) (found live: the poisoned chik2 run had 19/20 seeds
  cleanly shifted plus one at k ~ 1e153 and PASSED the z gate at max
  |z| = 3.27). Disjoint ranges are exact under exchangeability
  (P = 2/choose(2n, n) ~ 1.4e-11 at 20 seeds); the check self-disables
  below ~10v10 seeds, so quick runs are unaffected.

Fail-on-poisoned proofs (each: poison -> R CMD INSTALL --preclean into
the private .Rlib -> targeted compare vs a clean-build full local
baseline (/tmp/new-scen-full.rds, NOT recorded into benchmarks/) ->
revert -> reinstall):

Poison 16 (model.hpp:721 `noise = residualVariance / w` -> drops / w):
  wtgp 39 summaries, max |z| = 24.93, 29 with |z| > 3,
  26 with |z| > 4 <- FAIL (worst: vprop.2 z=24.93, k.mean z=-20.31,
  fhat.test.5 z=-18.52, ...) exit 1

Poison 12 (model.hpp:2334 `weightScratch[j] += w` -> `+= 1.0`):
  grouped 47 summaries, max |z| = 1266.61, 14 with |z| > 4 <- FAIL
  (worst: sigma.mean z=-1266.61, tau.mean z=-584.5, sigma.sd
  z=-158.41, tau.sd z=-110.59, vprop.2 z=51.75, ...) exit 1
  The mis-scaled group precision blows tau/ranef up; with the gamma
  prior the gate z-fails in ~2 s instead of hanging (the tinytest
  half-Cauchy hang the sweep recorded).

Poison 8 (model.hpp:1729 shape `0.5 (M + nu)` -> `0.5 (M + 2 nu - 1)`):
  chik2 39 summaries, max |z| = 3.27, 1 with |z| > 3,
  1 with disjoint seed ranges <- FAIL
    disjoint ranges: k.mean (4.14 vs 8.07e+151)  exit 1
  (Poisoned per-seed k.mean: 19 seeds at 5.2-6.3 vs baseline 3.9-4.4 -
  zero overlap - plus one diverged seed at 1.6e153 that degenerates the
  Welch z; the disjoint-range channel was added for exactly this.)

After the final revert + --preclean reinstall, the targeted compare of
wtgp/grouped/chik2 against the clean local baseline reports
"identical draws (same RNG stream)" for all three (exit 0).

### Sub-item 1: benchmarks/R/bcf-exact-weak.R (poison 13)

Implemented per the unblocked design (bcf-weak SBC: prior-dominant BCF
at n = 40 is well-behaved). K = 2 cells, n = 40, weak mu
(muTrue +-0.2, noise sd 1), a free with sd.control = 0.5, b fixed
(0, 1); exact side reuses bcf-exact.R's machinery with the a integral
under a = sd.control * tan(t) (Cauchy prior uniform in t, no tail
truncation); observables E[a mu], E[tau] (tolerance 0.02) and
P(|a| <= 0.25/0.5/1) (tolerance 0.025), the spread channel a mean-only
check is blind to. |a| lag-1 autocorrelation of the thinned draws is
reported per seed (clean: ~0.04 at thin 15, so the conservative
budget holds; under the poison it rises to 0.4-0.7). Full mode < 1 s.

Clean run (full): E[a mu] gap 0.0009, E[tau] gap 0.0010,
P(|a|<=q) gap 0.0041; exit 0.

Fail-on-poisoned proof (poison 13: chain.hpp:2049
`double aPrec = 1.0 / bcf_->aVariance` -> `double aPrec = 0.0`;
--preclean into the private .Rlib):
  |a| lag-1 autocorrelation of thinned draws, per seed: 0.41 0.73
  E[a mu]        exact -0.0556 -0.0400 | sampler -0.0884 -0.0849 (gap 0.0450)
  E[tau]         exact +0.1065 +0.1030 | sampler +0.1439 +0.1420 (gap 0.0390)
  P(|a|<=0.25,0.5,1) exact +0.4199 +0.7271 +0.9238 | sampler +0.0327
    +0.0650 +0.1030 (gap 0.8208 <- FAIL)
  FAIL: BCF a-glue posterior deviates from exact (exit 1)
Poison reverted (git checkout src/bartcore/chain.hpp), --preclean
reinstall, gate re-passes (quick, gaps 0.0019/0.0026/0.0166, exit 0).

## Status

- 2026-07-09: ALL SIX SUB-ITEMS DONE AND VERIFIED (each committed the
  moment its fail-on-poisoned proof and clean re-pass landed). Nothing
  half-built. Sub-item status:
  1. bcf-exact-weak.R - done + verified (poison 13 FAIL, clean pass).
  2. wtgp equivalence scenario - done + verified (poison 16 FAIL).
  3. grouped equivalence scenario - done + verified (poison 12 FAIL).
  4. linear-exact.R - done + verified (poison 15 FAIL, clean pass).
  5. bd-balance.R - done + verified (poison 2 FAIL, clean pass).
  6. chik2 equivalence scenario - done + verified (poison 8 FAIL via
     the new disjoint-range compare channel).
  Next concrete step: NONE for the implementer - ready for landing.
  The orchestrator re-records the extended 21-scenario equivalence
  anchor at landing (no baseline was recorded here by instruction).

## Final clean gates (worktree tree clean of engine edits)

- R CMD INSTALL --preclean into the private .Rlib: OK (the final
  install after the last poison revert).
- tinytest::test_package("dbarts"): 2498 passed, 0 failed.
- Equivalence compare vs benchmarks/baselines/equivalence-cf99a00.rds
  (run in two EQUIVALENCE_SCENARIOS halves to bound call length): all
  18 existing scenarios "identical draws (same RNG stream)"; wtgp,
  grouped, chik2 reported "not in baseline" with a warning; exit 0
  both halves.
- linear-exact.R full: all |z| <= 1.5, exit 0. bd-balance.R full: all
  |z| <= 1.3, exit 0. bcf-exact-weak.R full: gaps 0.0009/0.0010/0.0041,
  exit 0.
- air format + lintr clean on all touched R files (the single
  object_usage lint on equivalence.R's bartcore_bart predates this
  work, present at HEAD).
- tests/cpp untouched (no C++ test changes in this item); engine
  sources bit-identical to HEAD after every poison revert (verified
  via git status/diff each cycle).
- Private library note: all installs/runs used R_LIBS=<worktree>/.Rlib
  (~/.Renviron overrides R_LIBS_USER); the .Rlib directory is
  untracked scratch, not to be committed.
