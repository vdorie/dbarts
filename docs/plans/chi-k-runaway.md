# chi-k-runaway

agent: opus
rng: neutral - the cap is applied after the gamma draw, so RNG consumption
  is unchanged and every draw below the sentinel is bit-identical to the old
  code (verbatim FP order preserved). No shipped scenario engages the cap, so
  the equivalence suite is byte-identical.
budget: ~11 lines engine, one component test case, ~4 tinytest checks, one
  man-page sentence.

## Goal

Bound the sampled end-node precision k (src/bartcore/model.hpp,
ChiKHyperprior::draw) so that an improper or weak prior scale cannot let the
k Gibbs draw run away, while staying bit-identical for every run where the
bound never engages.

## Context

chi()'s default scale = Inf (R/model.R) is an improper prior on k. When
leaves are prior-dominated the k-Gibbs fixed-point growth factor
(1 + df / numLeaves) exceeds 1, so legal settings run away silently: chi(100,
Inf) does so deterministically (k -> 1e25+, leaf sd -> 0, the forest fits
nothing, sigma -> the marginal sd, no error raised); few-tree df = 1.5 does so
stochastically (k -> 1e21 observed). The finite-default-scale question is out
of scope (a separate research item); the decided remedy is a sentinel cap on
the sampled k alone, no warning path.

## Fix

A static constexpr ChiKHyperprior::maxDraw = 1e6. draw() computes the sqrt of
the gamma exactly as before, stores it, and returns min(k, maxDraw) via a
single `k > maxDraw ? maxDraw : k`. The gamma draw and its arithmetic are
untouched, so any run whose draw stays below 1e6 evaluates the identical FP
sequence and is bit-identical to the old code. At k = 1e6 the leaf sd is
already statistically zero on the standardized [-0.5, 0.5] response, so the
cap is behavior-neutral outside the runaway regime and bounds both runaway
modes (deterministic and stochastic).

State-restore decision (setState / installForest, src/bartcore/chain.hpp):
NOT capped there, deliberately. draw() is the sole stochastic origin of a
runaway k; with the cap in place no freshly run sampler can ever save a
runaway k, so a state produced by any supported workflow already carries a
k <= 1e6. setState only copies the saved k back (forest.k = fs.k), and a
pre-1.0 saved state is not a compatibility target (the 1.0-0 cutover reset the
compat contract). Capping in the two restore paths would add code to silently
rewrite a user's explicitly saved k for a case that no supported input can
produce, against the minimal-diff and cap-the-Gibbs-draw scope. So the cap
lives only where the runaway is generated.

## Test

- tests/cpp/test_model.cpp testChiKHyperprior gains two cases: (1) a
  prior-dominated draw (df = 1000, sumSq = 1e-30, leafScale = 1e-6, so the
  gamma mean ~ 1e21 and sqrt ~ 3e10) returns exactly ChiKHyperprior::maxDraw
  on 1000 consecutive draws; (2) a healthy draw (df 1.5, infinite scale, the
  case's leaf inputs) computed off a privately seeded rng equals the inlined
  sqrt-of-gamma formula bit for bit (checkNear tolerance 0.0) and sits below
  the cap. The private rng keeps the shared model stream untouched.
- inst/tinytest/test-binaryResponse-hyperprior.R gains a runaway block: a
  binary bart(k = chi(100, Inf), ntree = 5) fit runs to completion with every
  sampled k finite, all k <= 1e6, and at least one k == 1e6 (the cap engaged);
  plus a healthy bart(k = chi(1.5, Inf)) fit at the default tree count whose k
  trace exists and stays entirely below 1e6 (the cap is behavior-neutral).
  bart()'s own default is a fixed k = 2 (only bart2 defaults to the sampled
  chi(1.5, Inf)), so the healthy case passes the default chi hyperprior
  explicitly.

## Verification

- man/dbartsPriors.Rd chi() entry: one sentence noting the sampled k is capped
  at 1e6, the regime an improper scale can otherwise drift toward. ASCII-clean.
- Gates (worktree, private library
  /Users/vdorie/.claude/jobs/7fe13675/tmp/Rlib-chik):
  1. R CMD INSTALL --preclean -l <lib> - DONE (dbarts).
  2. tests/cpp clean rebuild, ./test_bartcore - exit 0, 93 ok lines, "all
     tests passed", including "ok: chi-k hyperprior".
  3. tinytest::test_package("dbarts") - 2530 passed / 0 failed (2526 baseline
     + 4 new checks).
  4. Equivalence vs benchmarks/baselines/equivalence-de67cbb.rds - all 21
     scenarios "identical draws (same RNG stream)", 21 compared / 0 skipped
     (both chik scenarios identical: the cap never engages in the suite).
  5. air format --check and lintr on the one touched R file
     (inst/tinytest/test-binaryResponse-hyperprior.R) - clean, 0 lints.

## Status

- 2026-07-10: fix + tests + doc landed on wt/chi-k-runaway; setState left
  uncapped per the recorded reasoning; all gates green.
- 2026-07-10: LANDED as 4797bc0 (squash of wt/chi-k-runaway).
  Reviewer approved the setState no-cap reasoning and re-ran the
  gates independently on the landed tree from a private library:
  preclean install clean, tests/cpp clean rebuild all pass,
  tinytest 2530/0, equivalence 21/21 identical draws vs
  equivalence-de67cbb.rds.
