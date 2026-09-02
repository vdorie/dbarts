# test-fit-parallel

agent: opus
rng: neutral (order-preserving by construction)
budget: ~100 lines

## Goal

Test-fit routing parallelizes over rows when numTest is large enough to
matter; results are byte-identical to the serial path at any thread
count.

## Context

- Serial today: findBottomNodeForRow per row per tree
  (src/bartcore/chain.hpp; [[tree.hpp:808@13e6154c]]). Rows are independent and
  writes are disjoint (each row owns its output slot), so partitioning
  rows across the existing chain worker threads preserves order and
  values exactly.
- Matters when numTest ~ n (TODO note); common in xbart folds and
  bart2 test-set workflows.

## Constraints

- Threads: reuse the chain-worker pool budget (min(numThreads, ...));
  no new threading machinery - a std::thread row-range split or the
  within-chain-threading pool if that lands first (note the
  interaction; do not block on it).
- Threshold: serial below a measured numTest cutoff so small test sets
  pay nothing.
- Out of scope: training-fit routing; predict()'s saved-tree replay
  (same pattern, separate call site - include only if trivially shared).

## Steps

1. Row-range split of the test-fit loop behind numThreads > 1 and
   numTest >= cutoff (measure the cutoff via benchmarks/).
2. Component test: multi-threaded test fits bitwise-equal serial ones
   (mirrors the existing chain thread-count invariance test).
3. Bench: record the win at numTest = n = 1e5.

## Verification

- Component tests; full tinytest; equivalence exact.
- bench-sampler compare: no regression at small numTest; recorded
  improvement at large.

## Landing note (2026-07-07, see merge commit)

Landed: routeTestRows splits test-fit routing across the chain's
share of the thread budget via a lazy per-chain misc_mt pool
(rebuilt on budget change, nothing allocated below the cutoff);
all six per-row routing loops converted; setNumThreads now
propagates to chains so setControl reaches the routing decision.
Cutoff 65536 (2-thread break-even ~65k); measured 1.48x wall at
numTest = n = 1e5 single-chain, 3.38x on the routing kernel at 8
threads / 262k rows. Byte-identical at any thread count (component
invariance test, 3e5 rows, 1/2/8 threads). Gates: full C++ suite +
fuzzer, tinytest 2470 ok, equivalence exact 18/18. bench-sampler
compare deferred: maintainer gate, expected no-op (its scenarios'
n.test = 25 sits far below the cutoff, leaving the serial path
untouched).
