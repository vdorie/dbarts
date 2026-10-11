# gp-sort-path: a gp leaf sorts its members through a bitmap, and the speed bound is read over more seeds

Status: PLANNED 2026-10-11, at bartcore f7fe0230. No open call; the first of the
[Calls made in planning](#calls-made-in-planning) is one the maintainer may want to mark.

agent: opus implementer, one; blind critique of this plan first; one opus reviewer who runs the mutants
below. Tier: changes draws ([Process by risk](README.md#process-by-risk)): the class is neutral, but the
sorted list feeds every gp draw, so a wrong sort moves them all.
rng: NEUTRAL, bit for bit, for every sampler. Only how
[`GPGaussianLeaf::sortedMembers`](../../src/bartcore/model.hpp) produces its array changes, never the array.
window: before the merge (dec-A192); engine slices stay serial. Either side of response-scale-rows's fix
round: as its branch stood on 2026-10-11 it shares model.hpp and test_model.cpp with this slice in disjoint
hunks (its model.hpp hunks all lie in the response models after `GPGaussianLeaf`), not test_peer.hpp.
budget: ~190 lines (model.hpp ~35, test_peer.hpp ~20, test_model.cpp ~115, records ~20).

## Goal

A gp leaf under the size cap gets its members in observation order in about half the time the sort takes
today, with every draw unchanged. The bound dec-A192 set (the member-order slice at or under 1.05 of its
base, warm, on every cap-64 arm) is read where the two builds' diverging chains cannot move it. No state
format, C API, facade or R surface change.

## Context

The item says cap64-n500-t50 reads 1.05 to 1.08 warm on the quiet x86-64 box (A = da19dec4, B' = 734441f1)
and asks for a cheaper sort. Measured with counters in scratch copies of both commits (scratch/gpsp:
replay.R, final.R, work.R) on a Mac under a load of 5 to 13: every count is exact, every time indicative.

- **Most of the 1.05 to 1.08 is not code.** The harness (scratch/gpimpl/ab/ab2.R) restores both builds
  from ten states a seed and times one cold sweep and two warm ones. The builds draw different values from
  the first gp draw on, so each of the thirty paths differs between A and B', and is the same pair of
  paths in every repetition and every rerun. On the report's three seeds B's warm sweeps hold more gp work
  than A's: members 1.025, their squares 1.040, their cubes (the Cholesky factors) 1.049. By seed the cubes
  read 1.019, 0.984 and 1.148, and the box's four runs 1.04 to 1.08, 1.01 to 1.03 and 1.10 to 1.12. Over
  24 seeds the cubes ratio has mean 0.994 and sd 0.062 and pools to 0.992: chance, not a property of
  either build; the gp path's time then reads B'/A 1.012 and 1.029 in two passes (seeds 1 to 3 alone 1.08
  and 1.09). These are the box's own paths: its csv files and the replay carry the same gp shares to
  fifteen digits, for both builds and all three seeds.
- **The sort is the one cost of any size the member-order slice added**: about 2.5 percent of a warm
  sweep on this arm (table). 96 percent of gp computations reach it. The spans are in random order (0.6
  descents an element; a random permutation has 0.5), and 99 percent of the spans the in-order test passes
  hold under 8 members. The draw cache's member copy, the lookup's re-check of the sorted members and
  eviction on bottom-ness and count cost under 6 us a sweep together; the base pays the last two as well.
- **The slice also saved work.** A rejected birth leaves the span permuted
  ([`birthOrDeathMove`](../../src/bartcore/moves.hpp)), which cost the base its cached kernel and factor;
  keyed on sorted members the entry survives (3.8 cholK factorizations a sweep against 4.6).
- **Why this arm.** Not leaf sizes (sorted spans average 28 to 32 members in every arm) and not cache
  misses (hit rates 58 to 74 percent; this arm 66). It has the most gp computations a sweep, 80 against 21
  and 26, so the sort's share is the largest, 2.5 percent against 1.8 and 1.2; and its three seeds drew
  the work excess, where n500-t10's read 1.006 and n2000-t50's 0.985.

Per warm sweep on the report's thirty states (us: in-process timers):

| | n500-t50 | n500-t10 | n2000-t50 |
|---|---|---|---|
| leaf evaluations; of them gp (under the cap) | 301; 80.0 | 68; 20.8 | 294; 26.1 |
| gp computations that sort; mean span | 76.5 (96%); 28.4 | 20.1 (97%); 31.6 | 24.7 (95%); 29.2 |
| sorted spans under 16, 16 to 31, 32 to 64 members | 34%, 23%, 43% | 29%, 24%, 47% | 26%, 35%, 40% |
| kernel-cache lookups (32 members up); hits; builds | 33.0; 21.8; 11.2 | 9.5; 6.5; 3.0 | 9.7; 5.6; 4.1 |
| the base's builds; cholK factorizations B', A | 11.5; 3.8, 4.6 | 3.5; 0.8, 1.2 | 4.5; 1.1, 1.3 |
| gp work B'/A: members, squares, cubes | 1.025, 1.040, 1.049 | 0.983, 0.992, 1.006 | 0.982, 0.986, 0.985 |
| us in the gp path, A and B' (whole sweep, about) | 502, 545 (670) | 171, 174 (250) | 183, 183 (470) |
| us sorting; copying members; lookup re-check; eviction | 16.5; 0.7; 2.5; 2.5 | 4.6; 0.2; 0.7; 0.6 | 5.5; 0.2; 1.0; 1.7 |

## Design

All in [`GPGaussianLeaf`](../../src/bartcore/model.hpp); prototype in scratch/gpsp/proto (proto.diff).

1. [`GPGaussianLeaf::sortedMembers`](../../src/bartcore/model.hpp) keeps its in-order return of the span
   itself, then fills `memberScratch_` by one of two methods, chosen by a private static predicate
   `sortsByBitmap(numObs, numWords)`:
   - the bitmap, when the span has at least 8 members and the bitmap at most 4 words a member: set each
     member's bit in a private `mutable std::vector<std::uint64_t> memberBits_` (one bit an observation),
     then walk the words in order, appending each set bit's index (`std::countr_zero`, then clear the
     lowest bit) and zeroing the word once read;
   - otherwise today's copy and `std::sort`.
2. `memberBits_` is sized `(numObservations_ + 63) / 64` zero words beside `u_` in
   [`GPGaussianLeaf::reinitialize`](../../src/bartcore/model.hpp), and is all zero between calls: every set
   bit lies in a word the walk visits and clears. Its threading contract is `memberScratch_`'s (the owning
   chain's thread; the test-row pool reads `memberBuffer_` and never sorts).
3. `friend struct TestPeer;` in the leaf. Callers, the cache key, eviction and the draw cache do not change.

Why bits cannot move: a leaf's span holds distinct observation indices below `numObservations_` (the
tree's spans partition the rows, and the leaf already indexes `u_` by them). Distinct integers have one
ascending arrangement, so the walk, `std::sort` and the base return the same array, and every caller reads
the array's contents only. The two thresholds choose a method, never a value.

Candidates, timed on 4096 random spans of 4 to 64 members (scratch/gpsp/sortbench.cpp; ns a span):

| candidate | effect | why bits cannot move, or why not taken |
|---|---|---|
| bitmap walk (taken) | n 500: 362 to 121; n 2000: 367 to 199; n 8000: 361 to 287; n 32000: 358 to 448, hence 4 words a member; spans of 4 to 8 at n 500: 37 to 50, hence 8 members. In place at n500-t50: 16.5 to 7.3 us a warm sweep | the one ascending arrangement |
| insertion sort under a size | 397 against 362; 687 against 623 at 48 to 64 members | same array; no gain on spans in random order (13 inversions an element) |
| cheaper in-order test | none: `std::is_sorted` leaves at the first descent, and 4 percent of spans are in order | same array |
| sort once a leaf a sweep (a memo by node, or the cache entry's list from 32 members) | at most the draw's sorts, 39 percent | exact only with a membership check, itself a pass over the members against a table of n (the walk's own work), or a hook in the moves; a lever on a miss |
| members kept in order at birth, change and swap (a stable partition) | removes the sort | NOT neutral: every leaf model sums in span order ([`Tree::partitionByPredicate`](../../src/bartcore/tree.hpp)), so every fit's draws move; another tier, and the maintainer's question |
| drop the lookup's re-check or the member copy | 2.5 and 0.7 us of 545 | the re-check is what makes a hit exact |

Chosen: the bitmap walk alone. It takes about 55 percent off the sort on this arm, in one function, on
a proof that needs no test of values. Prototype against the plain build: nine gp fits `identical()` in
draws, test fits, a copy's continuation and the serialized state (scratch/gpsp/bits.R: 200 to 5000 rows;
caps 64, 100 and 256; zero weights; probit on two chains; three designated columns). On the harness, chains
identical, n500-t50 reads warm 0.984 and cold 0.986 pooled over 3 seeds x 5 reps (medians 0.991, 0.996).

## Steps

1. model.hpp, ~35 lines: the predicate, its two constants, `memberBits_` and its sizing, the bitmap branch,
   the friend line; `sortedMembers`' comment says both methods return the one ascending arrangement.
2. [test_peer.hpp](../../tests/cpp/test_peer.hpp), ~20 lines: `TestPeer` gains three readers,
   `sortedMembers(leaf, tree, node)`, `sortsByBitmap(numObs, numWords)` and
   `holdsKernelEntry(leaf, tree, nodeIndex)` (the tree's entry at that node has members).
3. [test_model.cpp](../../tests/cpp/test_model.cpp), ~115 lines:
   - new `testGPLeafSortedMembers`: leaves over 64, 110, 500, 4097 and 20000 rows; spans of 2, 7, 8, 9, 31
     and 64 members in random order from a fixed generator, each holding rows 0 and n - 1; every result
     equals a `std::sort` of the span; two different spans in a row (a word left set would show); an
     in-order span returns the span's own pointer; the predicate's table at 7, 8 and 9 members and at 4
     words a member and one more; a leaf reinitialized from 100 rows to 300 sorts a span of rows above 200.
   - [`testGPLeafMemberOrder`](../../tests/cpp/test_model.cpp), at 48 members: after the score builds the
     entry and the span is reshuffled, `beginTreeDraw` leaves the entry held. This pins without a clock
     the mutant the member-order review caught by timing alone (eviction against the unsorted span).
4. Records, ~20 lines: a dated paragraph in [gp-leaves.md](../design/gp-leaves.md) under
   [Stage 4 addendum: members in observation order (2026-10-09)](../design/gp-leaves.md#stage-4-addendum-members-in-observation-order-2026-10-09)
   (the two methods; how the bound is read); TODO's item; this file's Landing note; the ledger entry. No
   NEWS and no help: gp leaves are new in 1.0-0 and nothing a user reads changes.

Reviewer's mutants, each to fail the named check: the word not zeroed after its walk (two spans in a row);
`& 63` or `* 64` wrong, or bits read from the top (every equality); the word count without its round-up
(row n - 1 at 110 and 500 rows, under the address sanitizer); `memberBits_` not resized (the reinitialized
leaf, same); the predicate's comparison flipped or its floor dropped (the table); the in-order return
dropped (the pointer); eviction against the unsorted span (the held entry).

## Verification

Against the slice's own library, independently of the implementer:

- `cd tests/cpp && make && ./test_bartcore`, plain and with `OPT="-O2 -g -fsanitize=address,undefined"`
  under `ASAN_OPTIONS=detect_container_overflow=0`.
- `R CMD INSTALL --preclean -l <lib> .`; the full tinytest suite, 0 failures, with
  [test-gp-leaves.R](../../inst/tinytest/test-gp-leaves.R) and
  [test-mutate-then-serialize.R](../../inst/tinytest/test-mutate-then-serialize.R) unchanged; those two
  also on the R-loaded sanitizer path (a new buffer indexed by observation).
- Reference build (`--preclean --configure-args=--enable-reference-build`), counted per scenario with no
  `max |z|` line ([Gate hygiene](README.md#gate-hygiene)): `equivalence.R compare
  benchmarks/baselines/equivalence-2b48939b.rds --bitwise --strict-coverage` 55 of 55, where
  ["GP leaves with a chi hyperprior"](../../benchmarks/R/equivalence.R) and
  ["GP leaves under NON-UNIT weights"](../../benchmarks/R/equivalence.R) are the two that build gp leaves
  (300 rows, cap 100: the bitmap from 8 members, `std::sort` below); `bcf-equivalence.R compare
  benchmarks/baselines/bcf-equivalence-1b7d730c.rds --bitwise` 15 of 15 and `multinomial-equivalence.R
  compare benchmarks/baselines/multinomial-equivalence-80b1c8d4.rds --bitwise` 11 of 11 (no gp leaf); the
  four `test-reproducibility-*.R` files unchanged (no gp leaf). Nothing is re-recorded.
- scratch/gpsp/bits.R on the base and on the slice: nine of nine `identical()`. Outside tests/cpp it is
  the only check that reaches `std::sort` on wide data (5000 rows).
- `R CMD check --as-cran` from a clean tarball; lintr, air, rc-codoc, win-drift, doc-freshness.

Timing, on the quiet box, orchestrator-run: shipped builds of A = da19dec4, B' = 734441f1 and B = the slice
tip, each `--preclean` into its own library; one thread, one chain, cap 64, as the report ran it.

1. `ab2.R <A> <B> out.csv cap64`: the report's seven arms, seeds 1 to 3, 5 reps. Bound: n500-t50 warm at
   or under 1.05; no other arm more than 0.01 above its figure in the report.
2. `ab2.R <B'> <B> out.csv cap64`: the states are burned in on B', so the chains are identical (the csv's
   two share columns agree in every row) and the ratio is the code alone. Bound: n500-t50 warm under 1.00
   (expected 0.975 to 0.99); no arm above 1.01, warm or cold.
3. scratch/gpsp/ab3.R (ab2.R with the seeds and reps as arguments) on `cap64-n500-t50`, 24 seeds, 2 reps:
   A against B, and A against B'. Bound: A against B at or under 1.05 warm. Expected from the Mac: B'
   1.01 to 1.02, B about 1.00.

Reading 1 can miss with the code sound: its three seeds carry a work ratio of 1.049 in the gp path, about
three quarters of the sweep, and the Mac's gain is 1 to 2.5 points. A miss of 1 with 2 and 3 passing meets
the bound; the Landing note records all three. If 2 misses, run scratch/gpsp/sortbench.cpp on the box and
retune the two constants. If 3 misses with 2 passing, stop and report: the excess is then not the sort.

## Budget and stops

Planned ~190 lines; recent engine slices ran about twice their plans, mostly in tests, so forecast 280 to
400. Implementer about two hours; gates about three hours of machine time; the box about forty minutes.

Stop and report when: the diff passes ~400 lines; an equivalence scenario or a snapshot moves, or bits.R
differs (the invariant is broken, not a baseline stale); a path is found that hands a gp leaf a repeated
member or one at or above `numObservations_`; the change needs anything outside `GPGaussianLeaf` and the
tests; reading 3 misses.

## Calls made in planning

- **Land the code although most of the miss is the harness.** The sort is 2.5 percent of a warm sweep on
  the densest gp arm, and the bitmap takes back 1 to 2.5 points for 35 lines and a proof. The alternative
  closes the item on reading 3 for B' alone, with no code. The cost of landing is one engine-tier round
  (critique, implementer, review with mutants, the gates above) for a gain of that size on gp fits only.
  Flagged because the item was written as a regression to fix and is mostly a measurement to correct.
- **dec-A192's bound is read over 24 seeds when the three-seed reading misses.** The alternative holds it
  to the report's three seeds, where the outcome turns on a work excess no code reaches (1.148 on seed 3).
  The cost is twenty minutes of box time and a bound whose reading changed after it was set; the
  three-seed figure stays in the record beside it.
- **The bitmap, not a memo or a stable partition** (the candidates table). The cost is one bit an
  observation per gp forest per chain (125 KB at a million rows) and a second sorting method to maintain.
- **Test hooks through `TestPeer`, not a counter in the leaf.** The alternative carries a kernel-build
  count beside the fallback tally, state no production path reads. The cost is a friend line.
- **bench-sampler.R compare is waived**, as dec-A192 waived it for the member-order slice: none of its
  arms builds a gp leaf. The cost is that the hot-path rule's own gate does not run; readings 1 to 3 do.
- **The two constants come from the Mac's microbenchmark**, not a tuning run on the box. The cost is a
  crossover that may sit elsewhere on x86-64; reading 2 shows whether it costs time.

## Not verified

- Nothing on x86-64: what `std::sort` and the bitmap walk cost there, and so the gain on the box. Every
  time in this file is the Mac's under load; the split of the box's 1.068 into work, sort and saved
  factorizations uses the box's work counts and the Mac's time shares.
- The prototype ran bits.R and the harness only: not tests/cpp, the tinytest suite, the equivalence trio,
  the snapshots or a sanitizer.
- That no path hands a gp leaf a repeated or out-of-range member: read from the partition code and the
  leaf's own indexing of `u_`, not tested.
- The planned tests are unwritten and their line counts estimates; response-scale-rows's fix round may
  add hunks; the blind critique has not run.
