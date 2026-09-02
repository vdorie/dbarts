This memo records the bench arm for TODO item monotone-leaf-branch-fill.
Outcome: not adopted (dropped); the variant is preserved on branch
archive/monotone-branch-fill at 98067fb3. See TODO for the current entry.

# monotone-leaf-branch-fill: bench report

Date 2026-08-26. Base tip de67cf07 (bartcore). The variant was developed on
branch `wt/monotone-bench` and is preserved verbatim on branch
`archive/monotone-branch-fill` at 98067fb3 (not adopted; see TODO).

Host: Apple Silicon, `sysctl -n hw.ncpu` = 10, R 4.6.1, clang -std=gnu++20 -O2.

## 1. Scope from the record

TODO entry `monotone-leaf-branch-fill` traces to the wave-5a fix-queue ledger
item 14 (`wave5-review-followups`, retired at 39fa379a), which cites the
wave-5a critique's item N6 (recorded in release-candidate-review.md's wave-5a
landing note): the rank is evaluated for EVERY leaf, so leafHasNoWeight's
first-vetoed-leaf early exit is gone, and the monotone rank scan is NEW work
([[moves.hpp:55-58@39fa379a]] returns first today) - a masked bench-sampler arm was called
for.

The empty-leaf veto slice added a per-leaf rank loop to
`logLikelihoodForBranch` (`[[src/bartcore/moves.hpp:67-96@39fa379a]]`) that runs for every
leaf model, monotone included. It fills the branch's bottom-leaf list into
`tree.bottomScratch`. Immediately afterwards, for a `ParamScoringLeafModel`,
`MonotoneConstantGaussianLeaf::logLikelihoodForBranchWithParams`
(`[[src/bartcore/model.hpp:729-771@39fa379a]]` at base) re-derived the SAME list into
`scratch.branch` with a second `tree.fillBottom(branchIndex, ...)`. That is the
"double branch-fill on the scoring path": one `fillBottom` per scored branch is
pure duplicate work, paid on both the current-branch and proposal-branch scores
of every birth/death of a monotone forest.

Nothing of this exists behind a flag; the duplicate fill is unconditional on
the monotone scoring path at de67cf07.

Invariant to preserve: bitwise-identical draws. The removed work touches no
arithmetic and no RNG stream, so the equivalence trio must report identical
draws, not merely statistical agreement.

## 2. The variant (diff summary)

3 files, +31 / -17 lines (`git diff --stat`):

- `src/bartcore/model.hpp` (+14 / -9): `ParamScoringLeafModel` gains a
  `const std::vector<std::int32_t>& branch` parameter;
  `logLikelihoodForBranchWithParams` takes the caller's list instead of calling
  `tree.fillBottom(branchIndex, scratch.branch)`, and uses it both for
  `numLeaves` and as the skip set of the per-leaf marginal loop. The now-dead
  `branch` member is dropped from `MonotoneNeighborScratch`.
- `src/bartcore/moves.hpp` (+4 / -1): `logLikelihoodForBranch` hands `bottoms`
  (already filled for the rank loop) down to the score.
- `tests/cpp/test_model.cpp` (+13 / -7): `ParamScoreMock` signature and the
  four direct call sites updated to pass an explicitly filled branch list.

The whole-tree fill `tree.fillBottom(0, scratch.allBottoms)` is retained - the
neighbor geometry genuinely needs it and the tree differs between the two
scores of a move, so it cannot be hoisted or cached.

## 3. Correctness gate (variant)

| Gate | Result |
| --- | --- |
| `tests/cpp`: `make && ./test_bartcore full` | all tests passed |
| tinytest `test_package("dbarts")` | 7468 results, 0 fails, 7468 passes |
| `equivalence.R compare equivalence-736bfb05.rds --strict-coverage` | exit 0; **43** "identical draws (same RNG stream)" lines; coverage 43 compared / 0 skipped; no z computed (identical), so max abs z = 0 |
| `bcf-equivalence.R compare bcf-equivalence-6e3b9fb8.rds --strict-coverage` | exit 0; **12** "identical (all N channels)" lines; every BCF channel bitwise identical |
| `multinomial-equivalence.R compare multinomial-equivalence-4d9a3337.rds --strict-coverage` | exit 0; **11** "identical (all 5 channels)" lines |

Bitwise across the trio, as a pure perf change requires. Timing proceeded.

## 4. A/B protocol

Same machine, same session, arms ALTERNATING base -> variant, 8 rounds
(the requested minimum is 5; 8 were run for a better noise estimate).
Two private libs built with `--preclean`: one from `git archive de67cf07`
(base), one from the variant worktree (variant). `benchmarks/kernels` holds no
binaries in either tree. `R_LIBS` set per arm on every call.

Each arm of each round runs, in order:

1. `Rscript benchmarks/R/bench-sampler.R` (print mode, full: 7 reps,
   500 samples) - the standard arms, none of which touch the monotone leaf.
2. A monotone-arm script written for this report (bench-sampler ships none).
   It mirrors bench-sampler's plain-run scenario shape (fixed seed/data,
   50-iteration warm-up, median of 5 reps of a 100-iteration `sampler$run`)
   but constructs the sampler with `dbarts(x, y, monotone = <named +1
   vector>)`, so every birth/death proposal goes through
   `MonotoneConstantGaussianLeaf`'s branch marginal. Four scenarios at n in
   {2000, 5000}, p = 10, 5 or 10 constrained columns, 75 or 200 trees;
   6-14 s of timed work per scenario, ~42 s per arm.

Nothing else heavy ran concurrently (only the negligible 2-minute gh poll).

Reported below: per-round values, the per-arm min over rounds (min-of-mins),
the variant/base ratio of those minima, and the within-arm spread
(max/min - 1) as the noise floor.

## 5. Results

Values are ms per iteration except where the metric column says otherwise.

### Monotone arm (where the variant applies)

| scenario | arm | r1 | r2 | r3 | r4 | r5 | r6 | r7 | r8 | min | spread |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| mono-n2000-p10-c5-t75 | base | 11.970 | 12.010 | 12.080 | 12.010 | 12.090 | 12.040 | 12.060 | 11.940 | **11.940** | 1.3% |
| mono-n2000-p10-c5-t75 | variant | 12.060 | 12.010 | 12.060 | 11.980 | 12.030 | 12.000 | 12.060 | 12.000 | **11.980** | 0.7% |
| | ratio | | | | | | | | | **1.0034** | +0.34% |
| mono-n5000-p10-c5-t75 | base | 13.170 | 13.120 | 13.060 | 13.400 | 13.100 | 13.070 | 13.420 | 13.350 | **13.060** | 2.8% |
| mono-n5000-p10-c5-t75 | variant | 13.180 | 13.080 | 12.970 | 13.130 | 13.020 | 12.990 | 13.090 | 13.130 | **12.970** | 1.6% |
| | ratio | | | | | | | | | **0.9931** | -0.69% |
| mono-n2000-p10-c10-t75 | base | 21.090 | 21.190 | 21.010 | 21.160 | 21.190 | 21.200 | 21.210 | 21.190 | **21.010** | 1.0% |
| mono-n2000-p10-c10-t75 | variant | 21.090 | 21.150 | 21.140 | 21.080 | 21.160 | 21.020 | 21.170 | 21.010 | **21.010** | 0.8% |
| | ratio | | | | | | | | | **1.0000** | +0.00% |
| mono-n2000-p10-c5-t200 | base | 28.330 | 28.620 | 28.390 | 28.370 | 28.360 | 28.360 | 28.310 | 28.480 | **28.310** | 1.1% |
| mono-n2000-p10-c5-t200 | variant | 28.380 | 28.570 | 28.390 | 28.400 | 28.370 | 28.350 | 28.310 | 28.380 | **28.310** | 0.9% |
| | ratio | | | | | | | | | **1.0000** | +0.00% |

### bench-sampler standard arms (where the variant does not apply)

| scenario (metric) | arm | r1 | r2 | r3 | r4 | r5 | r6 | r7 | r8 | min | spread |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| run-n1000-p10-t75 | base | 0.176 | 0.178 | 0.178 | 0.178 | 0.176 | 0.178 | 0.184 | 0.178 | **0.176** | 4.5% |
| run-n1000-p10-t75 | variant | 0.178 | 0.178 | 0.176 | 0.178 | 0.178 | 0.176 | 0.178 | 0.178 | **0.176** | 1.1% |
| | ratio | | | | | | | | | **1.0000** | +0.00% |
| run-n1000-p10-t200 | base | 0.468 | 0.474 | 0.474 | 0.474 | 0.474 | 0.484 | 0.476 | 0.476 | **0.468** | 3.4% |
| run-n1000-p10-t200 | variant | 0.472 | 0.474 | 0.474 | 0.472 | 0.474 | 0.476 | 0.472 | 0.474 | **0.472** | 0.8% |
| | ratio | | | | | | | | | **1.0085** | +0.85% |
| run-n10000-p10-t75 | base | 1.526 | 1.534 | 1.542 | 1.558 | 1.538 | 1.512 | 1.526 | 1.530 | **1.512** | 3.0% |
| run-n10000-p10-t75 | variant | 1.552 | 1.532 | 1.544 | 1.522 | 1.532 | 1.532 | 1.624 | 1.618 | **1.522** | 6.7% |
| | ratio | | | | | | | | | **1.0066** | +0.66% |
| run-binary-n1000-p10-t75 | base | 0.228 | 0.228 | 0.228 | 0.228 | 0.226 | 0.226 | 0.228 | 0.228 | **0.226** | 0.9% |
| run-binary-n1000-p10-t75 | variant | 0.228 | 0.228 | 0.228 | 0.228 | 0.228 | 0.226 | 0.228 | 0.228 | **0.226** | 0.9% |
| | ratio | | | | | | | | | **1.0000** | +0.00% |
| embedded-offset-run1-n1000-t75 (ms/gibbs step) | base | 0.260 | 0.260 | 0.260 | 0.264 | 0.256 | 0.256 | 0.256 | 0.260 | **0.256** | 3.1% |
| embedded-offset-run1-n1000-t75 | variant | 0.260 | 0.256 | 0.260 | 0.256 | 0.260 | 0.256 | 0.256 | 0.260 | **0.256** | 1.6% |
| | ratio | | | | | | | | | **1.0000** | +0.00% |
| setPredictor-accept-n1000-t75 (ms/update) | base | 0.259 | 0.262 | 0.260 | 0.266 | 0.259 | 0.260 | 0.257 | 0.258 | **0.257** | 3.5% |
| setPredictor-accept-n1000-t75 | variant | 0.258 | 0.263 | 0.263 | 0.261 | 0.261 | 0.260 | 0.259 | 0.260 | **0.258** | 1.9% |
| | ratio | | | | | | | | | **1.0039** | +0.39% |
| setPredictor-reject-n1000-t75 (ms/update) | base | 0.148 | 0.147 | 0.148 | 0.145 | 0.148 | 0.147 | 0.150 | 0.150 | **0.145** | 3.4% |
| setPredictor-reject-n1000-t75 | variant | 0.151 | 0.147 | 0.147 | 0.147 | 0.149 | 0.150 | 0.151 | 0.146 | **0.146** | 3.4% |
| | ratio | | | | | | | | | **1.0069** | +0.69% |

Every ratio lies in [0.993, 1.009]; every within-arm spread is 0.7-6.7%. No
metric's variant/base difference exceeds its own arm's round-to-round spread.

### Load stamps (1/5/15-minute loadavg before and after each arm)

| round | arm | before | after |
| --- | --- | --- | --- |
| 1 | base | 5.00 12.97 10.30 | 3.75 11.26 9.83 |
| 1 | variant | 3.75 11.26 9.83 | 2.76 9.74 9.36 |
| 2 | base | 2.76 9.74 9.36 | 3.35 8.74 9.01 |
| 2 | variant | 3.35 8.74 9.01 | 3.16 7.80 8.64 |
| 3 | base | 3.16 7.80 8.64 | 2.60 6.88 8.24 |
| 3 | variant | 2.60 6.88 8.24 | 2.25 6.06 7.85 |
| 4 | base | 2.25 6.06 7.85 | 2.60 5.51 7.53 |
| 4 | variant | 2.60 5.51 7.53 | 2.59 5.02 7.22 |
| 5 | base | 2.59 5.02 7.22 | 4.13 5.09 7.11 |
| 5 | variant | 4.13 5.09 7.11 | 2.63 4.51 6.77 |
| 6 | base | 2.23 4.03 6.43 | 4.43 4.29 6.37 |
| 6 | variant | 4.43 4.29 6.37 | 4.21 4.27 6.23 |
| 7 | base | 4.21 4.27 6.23 | 5.40 4.66 6.26 |
| 7 | variant | 5.40 4.66 6.26 | 3.66 4.30 6.02 |
| 8 | base | 3.66 4.30 6.02 | 4.54 4.38 5.94 |
| 8 | variant | 4.54 4.38 5.94 | 3.77 4.18 5.77 |

**FLAG: every round exceeded the loadavg <= 2 abort threshold.** The host's
loadavg floor with the owner's desktop session up (Firefox, Docker Desktop,
WindowServer, several claude processes) never fell below ~2.2 during the run;
it was 2.07/2.97/5.57 before any benchmark work started. macOS loadavg counts
blocked as well as runnable threads and overstates contention here: at the end
of the run `top -l 2 -n 0` reported **82.6% idle** (11.9% user, 5.5% sys) on 10
cores, i.e. roughly one busy core - the benchmark's own single thread. Rather
than abort every round, the alternating A-B-A-B ordering was relied on to share
any drift between arms, and the within-arm spread (0.7-6.7%) is reported as the
empirical noise floor. The two arms' loadavg stamps interleave with no
systematic difference. Rounds 7 and 8 carry the run's only visible disturbance
(`run-n10000-p10-t75` variant at 1.624/1.618 against its own 1.522 minimum);
min-of-mins absorbs it.

### Machine representativeness vs. the canonical baseline

`bench-sampler.R compare benchmarks/baselines/bench-sampler-ab1dc52.csv` on the
BASE lib, and the same comparison using this run's base min-of-mins:

| scenario | canonical ab1dc52 | base min-of-mins here | ratio |
| --- | --- | --- | --- |
| run-n1000-p10-t75 | 0.174 | 0.176 | 1.011 |
| run-n1000-p10-t200 | 0.464 | 0.468 | 1.009 |
| run-n10000-p10-t75 | 1.524 | 1.512 | 0.992 |
| run-binary-n1000-p10-t75 | 0.222 | 0.226 | 1.018 |
| embedded-offset-run1-n1000-t75 | 0.248 | 0.256 | 1.032 |
| setPredictor-accept-n1000-t75 | 0.247 | 0.257 | 1.040 |
| setPredictor-reject-n1000-t75 | 0.136 | 0.145 | 1.066 |

The five run/gibbs arms land within 1-3% of the canonical baseline, so the
machine is representative. The two setPredictor arms - the smallest metrics in
the grid, 0.14-0.26 ms - sit 4-7% high, and a single `compare` invocation
therefore exits 1 with two REGRESSION flags. That offset is present in BOTH
arms of this A/B (base and variant min-of-mins differ by 0.4-0.7% there), so it
is a host-level level shift, not a signal about the variant. It is NOT a claim
that tip de67cf07 regressed those arms; it is what an ab1dc52 comparison looks
like on this host in this session.

## 6. Why the effect is this small (profile)

`sample <pid> 12` on the base lib during `mono-n2000-p10-c5-t75`
(9645 thread samples):

- `metropolisJumpForTree<MonotoneConstantGaussianLeaf>` 5614 samples, of which
  `logLikelihoodForBranch` 3232 (58%);
- inside that, `twoLeafCoupledLogMarginal` 1675 + 1552 = 3227 - essentially all
  of the branch score - and inside that, `monotoneIntegrate` /
  `monotoneAdaptiveSimpson` over `coneProbability`, bottoming out in
  `Rf_pnorm5`. The adaptive-Simpson recursion is the monotone scoring path's
  cost.
- `Tree::fillBottom` appears out-of-line in the profile and accounts for
  **17 of 9645 samples = 0.18% of the thread**, summed over ALL of its call
  sites (moves.hpp's rank fill, the model's branch fill, the model's
  whole-tree fill, the leaf-draw fills in chain.hpp).

The variant removes one of those call sites, so the theoretical ceiling on the
win is a fraction of 0.18% - an order of magnitude below the 0.7-2.8%
round-to-round noise. The measurement and the profile agree.

## 7. Recommendation

**Do not adopt on perf grounds; inconclusive-by-construction as a speedup.**
The measured effect is 0.00% to -0.69% on the four monotone scenarios
(min-of-mins ratios 1.0034, 0.9931, 1.0000, 1.0000) and 0.00% to +0.85% on the
seven bench-sampler arms the change cannot touch - the two ranges overlap
completely, which is the signature of no effect. The round-to-round spread
within a single arm is 0.7-2.8% on the monotone scenarios and up to 6.7% on the
standard arms, so nothing here clears noise in either direction, and the
profile explains why: the eliminated `fillBottom` is at most 0.18% of the
thread even counting every other call site, against a scoring path dominated by
adaptive-Simpson quadrature over `Rf_pnorm5`. N6's cost concern is real in
structure but negligible in magnitude. The diff is nonetheless bitwise-clean on
the full trio and passes both engine gates, and it removes a genuine
duplicate computation and one dead scratch vector, so if VD wants it it should
be taken as a readability/tidiness change with a signature widening
(`ParamScoringLeafModel` grows a parameter, touching `tests/cpp/test_model.cpp`)
- not as perf. If the monotone leaf ever needs to get faster, the target named
by this profile is `coneProbability`'s quadrature (3227 of 3232 branch-score
samples), not the fills.

## Artifacts

- variant: preserved on branch archive/monotone-branch-fill at 98067fb3
  (not adopted; see TODO)
- supporting scripts, per-round logs, and the profile capture were left in
  the benching worktree's local scratch area; not part of the tracked record
