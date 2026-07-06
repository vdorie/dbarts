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
  (src/bartcore/chain.hpp; tree.hpp:808). Rows are independent and
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
