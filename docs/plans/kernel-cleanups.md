# kernel-cleanups

agent: sonnet
rng: neutral
budget: ~200 lines

## Goal

The fast serial moment accumulators are callable without the null-
thread-manager indirection, and misc.a no longer calls into R directly
(injectable output hook), closing the two cleanup candidates
kernel-vocabulary.md records.

## Context

- bartcore calls misc_htm_compute* with (nullptr, 0) everywhere
  (src/bartcore/tree.hpp:471-513, chain.hpp:548, model.hpp:1622-1624):
  the htm entry points with a null manager ARE the fast path
  (docs/design/kernel-vocabulary.md:92-98); the plain misc_compute*
  variants are several times slower.
- misc.a prints via Rprintf (htm), which benchmarks stub; the doc
  flags an injectable output hook as the fix.
- The thread managers themselves stay dormant (within-chain-threading).

## Constraints

- No behavior change; pure interface exposure. The htm entry points
  keep working (within-chain-threading revives them).
- constant-leaf-suffstat may replace some of these call sites; if both
  are in flight, land this first (it is mechanical) or rebase.
- Out of scope: deleting the slow plain variants' external symbols if
  anything links them (check benchmarks/kernels and tests/cpp).

## Steps

1. Expose the serial accumulator bodies as misc_compute*Fast (or
   promote them to the plain names if the slow variants prove
   caller-free) in stats.h; htm null-manager path delegates to them.
2. Switch bartcore call sites to the direct entry points; drop the
   (nullptr, 0) arguments.
3. misc output hook: a settable function pointer defaulting to Rprintf
   in the package build and fprintf(stderr) standalone; benchmarks/
   tests drop their stubs.

## Verification

- Component tests bitwise-identical results pre/post (same code paths,
  new names).
- Full tinytest; equivalence exact.
- benchmarks/kernels builds without the Rprintf stub.
