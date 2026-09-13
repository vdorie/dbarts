# valgrind-xbart

agent: opus (the measurement and the diagnosis).
rng: neutral. Nothing on a draw path is touched. The one change is to a
  test file's comparison operator, and only in the environment named
  below; the three bitwise equivalence baselines are untouched and expect
  IDENTICAL as before.
window: the pre-release sequence's memory check. Independent of the other
  pre-release arcs.
budget: ~25 lines in one tinytest file, plus this note.

## Goal

The nightly valgrind job is green at the tip, and the seven assertions
that failed under it on its last recorded run are explained: six in
["expectSameSweep"](../../inst/tinytest/test-xbart-reproducibility.R)
and one convergence-diagnostic snapshot in
["splitRhat"](../../inst/tinytest/test-convergence-diagnostics.R).
Neither is a defect in the package. Both come from one property of the
tool: valgrind emulates the x87 unit at 64 bits, so anything R accumulates
in a long double rounds differently under it than on the hardware.

## What valgrind reports at the tip

Nothing. The whole tinytest suite ran under memcheck with leak checking
full and origin tracking on, the way the nightly job runs it, against a
tip installed from a clean preclean build on an x86 box:

- definitely lost: 0 bytes in 0 blocks.
- indirectly lost: 0 bytes in 0 blocks.
- possibly lost: 0 bytes in 0 blocks.
- invalid reads, invalid writes, invalid or mismatched frees, conditional
  jumps on uninitialized values: none. The error summary is zero errors
  from zero contexts, and the log carries no instance of any of those
  words.
- still reachable at exit: about 305 MB in 64 thousand blocks, over 7.3
  million allocations and 17 GB allocated in total. That is R's own
  arenas and the session's live objects at quit, not a leak; the job's
  own criterion folds definitely- and indirectly-lost blocks into the
  error count and ignores this line.

8653 assertions ran. 8647 passed. The six that did not are the subject of
the next section, and they are the same six the nightly job reported. The
convergence-diagnostic snapshot that also failed there now passes, for
the reason given below.

One caveat on the coverage: the nightly job's container carries an R
built with valgrind instrumentation, which teaches memcheck about R's own
allocator and so reports uninitialized reads that a stock R would hide
behind a reused block. The box used here has a stock R, so this run is
the weaker of the two on that one axis. It is not the weaker one on
leaks or on invalid access, and the instrumented run reported the same
zeros.

## The six xbart assertions

The file asserts that a seeded `xbart` sweep returns at two, three and
four threads exactly what it returns at one. Six of those comparisons
failed under valgrind and nothing else in the file did.

The six are the one-worker-against-many ones, and only those. At one
thread [`xbart`](../../R/xbart.R) runs the (replication, fold) units in
the calling process; above one it runs them in `parallel` worker
processes. Valgrind traces the process it starts and not the children it
spawns, so a one-against-many comparison under valgrind compares an
emulated result against a native one. The comparison in the same loop
that does not cross that line - a rerun at a given thread count against
itself - passed, at every count, and so did the two that check the
caller's random stream is where it was.

Measured on an x86 box with the same source, the same options and the
same test file:

- Natively, every arm is bitwise identical: one thread, two threads and
  four threads return the same bits, and the file passes.
- Under valgrind, the two- and four-thread arms are still identical to
  each other and to all three native arms, bit for bit. The one-thread
  arm - the only one valgrind actually executed - differs from them in
  five of the eight reported cells, by one or two units in the last
  place: about 4e-16 on a loss near 2.4, a relative difference of 2e-16.

The sampler is not involved. A seeded fit's sigma draws, all six hundred
of its test predictions and all three thousand of its training
predictions are bit-identical under valgrind and natively, and the
runtime SIMD level the shared object dispatches to is the same (avx2) in
both, so no kernel and no reduction order inside the engine changed.
What changed is above it: `xbart`'s rmse loss is `sqrt(mean(...))` over
`rowMeans(...)`, and both of those R primitives accumulate in a long
double on x86_64. Probed directly in the same two processes, `mean()` on
a thousand normal draws differs by one unit in the last place between
them and `sum()/n` by sixteen.

So the property the test is defending - that no draw depends on which
worker ran a unit, or on how many there were - is intact. What is not
available under valgrind is the bitwise form of the comparison, because
the two sides are not running the same arithmetic.

## The convergence-diagnostic snapshot

One split-Rhat literal differed in the fourth decimal: 0.9567 expected
against 0.9574 observed, a relative difference near 8e-4, far too large
to be rounding and exactly the size the file's own comment predicts for
this fixture.

It is already fixed at the tip, by a commit that landed about an hour
after that run started. The rank-normalized Rhat folds the pooled draws
around their median; when the pooled count is even that median is the
mean of the two middle order statistics, and the two folded values built
from them tie mathematically. Whether the tie survives as a
floating-point equality depends on the width of the accumulator `mean()`
uses, and `rank()`'s average-ties rule then assigns a tied and an untied
pair different ranks, which moves Rhat at the 1e-3 level. The fix already
in the tree rounds every pinned fixture onto a 1/64 grid, so the two
middle order statistics sum exactly and no accumulator width can break
the tie. It was written for the difference between x86_64's long double
and arm64's plain double; valgrind is the same difference by another
route, and the guard covers it.

Confirmed by measurement: at the tip that file passes under valgrind, all
sixty-five assertions, and every one of its pinned diagnostics is either
bit-identical to its native value or within a few units in the last
place - nowhere near the 1e-12 the assertions allow.

## The fix

A test change, in
["expectSameSweep"](../../inst/tinytest/test-xbart-reproducibility.R)
only. The two one-against-many comparisons go through a small helper that
compares bitwise everywhere and to a 1e-12 tolerance when the process is
running under valgrind, detected by the preload valgrind puts in the
environment of the process it traces. Every other assertion in the file,
including the same-thread-count reruns and the caller's-stream ones,
stays exactly as it was.

The tolerance does not weaken the gate. A real thread-count dependence -
a unit seeded from its worker rather than from its own index, a fold
assignment that reads the worker count - moves a whole fit, which shows up
as a difference in the first significant digits of a loss, not in its last
bit. Twelve orders of magnitude separate the two.

Self-detection rather than a flag the workflow sets, because the tests
also run under `R CMD check`, and CRAN's own valgrind check would meet
this the same way with no workflow to set anything.

## Verification

On an x86 machine with valgrind, against an installed tip:

    R -d valgrind --vanilla -e 'tinytest::test_package("dbarts")'

expects every assertion to pass and a zero error summary with nothing
definitely or indirectly lost. The reproducibility file on its own, run
natively,

    tinytest::run_test_file("inst/tinytest/test-xbart-reproducibility.R")

expects twenty-four passing assertions, the one-against-many comparisons
still bitwise; run under valgrind it expects the same twenty-four, those
two to a tolerance.

## Status

Status: RUN AND REPORTED, 2026-09-13. The tip is clean under valgrind;
the xbart test change lands with this note; the snapshot needed no change.
