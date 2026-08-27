# tau-slice-stepout-cap

agent: fable
rng: neutral (integer guard short-circuits ahead of the untouched FP
  conditions; no RNG consumed unless the cap engages, and the cap never
  engages inside the gate suite)
budget: ~15 lines engine, one component test

## Goal

Bound sliceSampleOnce's stepping-out loops (src/bartcore/model.hpp) so
the grouped tau draw cannot hang when tau random-walks far from its
bracket, while staying bit-identical for every run where the bound
never engages.

## Context

Review-6 numerical finding (seed 4), REAL/high: only the shrinkage loop
had an iteration cap (1000); BOTH step-out while loops were unbounded.
Step-outs cost ~tau/width density evaluations: >1e5 iterations at
tau > ~2.5e5, an indefinite hang at tau > ~1e8. Reachable with legal
inputs: empty groups draw effects ~N(0, tau^2) with no mean-reversion,
so mostly-empty groupings let tau random-walk upward; heavy-tail
excursions at small J. Also the cause of the review-3 poison hang and
the review-4 half-Cauchy SBC intractability.

## Fix

Neal (2003)'s m cap on the step-out, applied per side: at most 10^4
width-expansions in each direction. Constant justification: healthy
runs measure under 100 expansions, so 10^4 is two orders of magnitude
of headroom while capping the worst case at ~2e4 density evaluations
(microseconds) instead of unbounded. When the cap engages the bracket
is simply the capped interval; it still contains the current point, so
the shrinkage loop samples correctly inside it. The integer counter
check (steps-- > 0) is prepended to each loop's && chain, so it
short-circuits before the existing floating-point comparisons: any run
where the cap does not engage evaluates the identical FP sequence and
is bit-identical to the old code.

## Test

tests/cpp/test_model.cpp testGroupedMath gains a stall-regime case: the
tau posterior at J = 5, sum b^2 = 5e24 (mode ~ sqrt(SS/(J+2)) ~ 9e11)
sliced from x = 0.45 with width 0.55 on (0, Inf). Unfixed, the right
step-out runs ~1.7e12 iterations - an indefinite hang. Fixed, it must
return a finite positive draw inside the capped bracket
x + (m + 1) * width. A dedicated seeded rng keeps the shared test
stream untouched.

## Verification

- Hang-on-unfixed: with the new test in place and src/bartcore/model.hpp
  reverted to the pre-fix version, `./test_bartcore model` (pty-recorded
  so output flushes per line) froze after "ok: linear leaf views" - the
  test immediately before grouped math - and spun at 100% CPU for over
  2.5 minutes before being killed. The fixed header completes the entire
  model suite in 1.21 s wall, so the stall case alone overshot the
  healthy suite by >120x and was still climbing toward its predicted
  ~1.7e12 iterations. Header restored, suite rebuilt and green.
- Gates (worktree, shared default R library):
  1. R CMD INSTALL --preclean . - PASS (DONE (dbarts)).
  2. tests/cpp clean rebuild, ./test_bartcore - 93 ok lines, "all tests
     passed", including the new case ("ok: grouped math").
  3. tinytest::test_package("dbarts") - 2498 passed, 0 failed.
  4. Equivalence vs benchmarks/baselines/equivalence-0e9ccca.rds - all
     21 scenarios "identical draws (same RNG stream)", 21 compared /
     0 skipped (cap never engages in the suite).
  5. No R files touched (diff is model.hpp, test_model.cpp, this doc).

## Status

- 2026-07-10: fix + test landed on wt/tau-slice-stepout-cap; hang
  verification recorded; all gates green.
