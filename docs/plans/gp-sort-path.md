# gp-sort-path: a gp leaf sorts its members through a bitmap, and its speed is read on identical chains

Status: PLANNED 2026-10-11 (built on the outcome of a rule fixed before its data,
[The rule and its outcome](#the-rule-and-its-outcome)). Revised to the blind critique of its first version,
whose nine findings are answered here, and to the x86-64 box's measurements. No open call.

agent: opus implementer, one; one opus reviewer who runs the mutants below. Tier: changes draws
([Process by risk](README.md#process-by-risk)): the class is neutral, but the sorted list feeds every gp
draw, so a wrong sort moves them all.
rng: NEUTRAL, bit for bit, for every sampler. Only how
[`GPGaussianLeaf::sortedMembers`](../../src/bartcore/model.hpp) produces its array changes, never the array.
window: before the merge; base d5b637ee or later. Engine slices stay serial; queued beside this one are the
variance-forest fix in [`Chain::setModel`](../../src/bartcore/chain.hpp), the forest slices and the
collapsed scale move. This slice touches `GPGaussianLeaf`, test_peer.hpp and test_model.cpp only; against
d5b637ee, response-scale-rows's model.hpp hunks all lie after `GPGaussianLeaf`, its test_model.cpp hunks
are its own tests and two lines of `runModelTests`, and test_peer.hpp is as it was.
budget: ~245 lines (model.hpp ~50, test_peer.hpp ~30, test_model.cpp ~140, records ~25).

## Goal

A gp leaf under the size cap gets its members in observation order for a fraction of what the sort costs,
with every draw unchanged: on the box, a warm sweep of cap64-n500-t50 at or under 0.98 of its base on
identical chains (the prototype read 0.941). No state format, C API, facade or R surface change.

## Context

The bound this item was opened on is not dec-A192's, which states none. It is the member-order plan's stop
condition, any A/B arm above 1.05 after the levers
([Stop conditions](gp-copy-continuation.md#stop-conditions); its accept was every arm at 1.03,
[Speed A/B](gp-copy-continuation.md#speed-ab)), restated as the rerun item's trigger, above 1.05 in any arm
([TODO:288-293](https://github.com/vdorie/dbarts/blob/7f948d0233606a6a0990e0e626f966228cccc3f9/TODO#L288-L293)).
The statistic is the pooled warm ratio: an arm's warm seconds summed over its seed x rep rows, B over A.

Box (quiet x86-64 Linux, g++ 13.3, shipped builds, one thread; scratch/x86speed/gpsp), cap64-n500-t50.
A = da19dec4, B' = 734441f1 (the member-order slice), P = f7fe0230, Q = P with the prototype:

| reading | seeds x reps | warm, pooled (per seed) | cold | at equal work: warm; cold |
|---|---|---|---|---|
| B'/A, the report's four runs | 1 to 3 x 5 | 1.050 to 1.068 | 1.021 to 1.035 | work itself 1.049; 0.999 |
| B'/A | 1 to 24 x 5 | 1.028 (0.949 to 1.097, se 0.008) | 1.021 | 1.033 (se 0.004); 1.019 (0.003) |
| B'/A, seeds not looked at before | 25 to 48 x 3 | 1.017 (0.916 to 1.172, se 0.013) | 1.015 | 1.021 (0.004); 1.015 (0.004) |
| Q/P, identical chains | 25 to 48 x 3 | 0.941 (0.914 to 0.964, se 0.003) | 0.934 | the same work |

Work is the sum over gp computations of the members cubed (the Cholesky factors), counted on the Mac in
instrumented copies of A and B' (scratch/gpsp: replay.R, box.R); the counts are the box's own paths, the
replay carrying the box csv's gp shares to 5e-16 for all 48 seeds and both builds. "At equal work" is the
per-seed regression of the log time ratio on the log work ratio, read at a work ratio of 1 (R2 0.76 and
0.90 warm).

- **What the three-seed miss was.** The harness restores both builds from ten states a seed and times one
  cold sweep and two warm ones; the builds draw different values from the first gp draw on, so each path
  differs between A and B' and is the same pair of paths in every repetition and rerun. Of the six warm
  points, two to three and a half are code (1.021 and 1.033 at equal work) and the rest is the three
  seeds' draw of paths (warm work 1.049; by seed 1.019, 0.984 and 1.148, where the box read 1.04 to 1.08,
  1.01 to 1.03 and 1.10 to 1.12). The cold column was code throughout: its work ratio is 0.999, and it reads
  1.02 to 1.035 on three seeds against 1.015 to 1.021 on 24. Seeds 1 to 24 had been looked at and pool to a
  work ratio of 0.992; seeds 25 to 48 pool to 0.994 and are the ones every reading below uses.
- **The code is the sort, less what the slice saved.** On the box the bitmap alone takes 5.9 points off a
  warm sweep and 6.6 off a cold one, so `std::sort` costs at least that there (the Mac's libc++: 2.5), and
  the member-order slice's kept kernels pay back three to four: a rejected birth leaves the span permuted
  ([`birthOrDeathMove`](../../src/bartcore/moves.hpp)), which cost the base its cached kernel and factor
  (3.8 cholK factorizations a warm sweep against 4.6).
- **Where the sort runs** (counts, the report's thirty states, per warm sweep): 80 gp computations of 301
  leaf evaluations at n500-t50 (21 of 68 at n500-t10, 26 of 294 at n2000-t50); 96 percent of them sort, the
  spans being in random order (0.6 descents an element; a random permutation has 0.5); sorted spans average
  28 members: 34 percent under 16, 23 percent 16 to 31, 43 percent 32 to 64. The draw cache's member copy,
  the lookup's re-check and eviction cost under 6 us of a 670 us sweep on the Mac, the last two in the base
  as well. This arm has the most gp computations a sweep; leaf sizes and cache hit rates match the others'.

## The rule and its outcome

Fixed on 2026-10-11 after the first version of this plan and its critique and before the box runs below
returned (scratch/gpsp-protocol.md): build if Q/P pooled warm on identical chains, seeds 25 to 48, is at or
under 0.98, or if B'/A at equal work on those seeds is above 1.04; otherwise record the figures and bring
"close with no code" to the maintainer. Outcome: Q/P 0.941, every one of the 24 seeds under 0.98; B'/A at
equal work 1.021. BUILD, on the first clause. The protocol's second Q/P run, with the libraries swapped,
had not been run when this was written. The three-seed figure stays in the record as history.

## Design

All in [`GPGaussianLeaf`](../../src/bartcore/model.hpp). The prototype (scratch/gpsp/proto.diff) is what
the box timed; it carries an earlier rule (8 members, 4 words a member) and one function.

1. Two private helpers, each filling `out` with `numObs` members in ascending order:
   - `sortMembersByBitmap`: set each member's bit in a private `mutable std::vector<std::uint64_t>
     memberBits_` (one bit an observation), then walk the words in order, appending each set bit's index
     (`std::countr_zero`, then clear the lowest bit) and zeroing the word once read. It asserts that
     `memberBits_` holds this data's word count and that the walk ended at `out + numObs`: a repeated member
     would otherwise leave a stale tail (tests/cpp keeps asserts; the shipped build defines NDEBUG).
   - `sortMembersByComparison`: today's copy and `std::sort`.
2. A private static predicate `sortsByBitmap(numObs, numWords)`: at least 16 members and at most 6 words a
   member. [`GPGaussianLeaf::sortedMembers`](../../src/bartcore/model.hpp) keeps its in-order return of the
   span itself, then calls the helper the predicate names, with `numWords` this data's
   `(numObservations_ + 63) / 64`.
3. [`GPGaussianLeaf::reinitialize`](../../src/bartcore/model.hpp) sizes `memberBits_` to that many zero
   words when the predicate holds at `min(maxLeafSize_, numObservations_)` members, and leaves it empty
   otherwise: the predicate grows with the member count, so no span under the cap can then pass (cap 64:
   above 24576 rows; cap 256: above 98304). It also clears `memberScratch_`. The bitmap is all zero between
   calls (every set bit lies in a word the walk visits and clears); its threading contract is
   `memberScratch_`'s.
4. `friend struct TestPeer;` in the leaf. Callers, the cache key, eviction and the draw cache do not change.

Why bits cannot move: a leaf's span holds distinct observation indices below `numObservations_` (the
tree's spans partition the rows, and the leaf already indexes `u_` by them). Distinct integers have one
ascending arrangement, so both helpers and the base return the same array, and every caller reads the
array's contents only. The predicate chooses a method, never a value.

The predicate, from the box's sortbench (scratch/gpsp/sortbench.cpp; ns a span: `std::sort`, insertion
sort, bitmap walk):

| rows (words) | 8 to 15 members | 16 to 31 | 32 to 64 |
|---|---|---|---|
| 500 (8) | 83.5, 82.6, 63.0 | 273, 202, 80 | 699, 536, 107 |
| 2000 (32) | 82.5, 82.4, 101 | 275, 202, 162 | 701, 535, 214 |
| 8000 (125) | 82.5, 82.1, 144 | 275, 203, 215 | 700, 534, 358 |

- Floor 16: from 16 members the walk beats `std::sort` in all six cells, by 22 to 85 percent. At 8 to 15 it
  wins 20 ns at 500 rows and loses 19 and 61 at 2000 and 8000. The floor gives up that 20 ns on about 14
  spans a sweep of n500-t50 (0.3 us of 900) and removes the losing low edge the critique found.
- 6 words a member: on the 8000 row the walk costs about 6 ns a member over 144 at 11.5, so about 170 at 16
  members, where `std::sort` interpolates to about 155; between the 2000 and 8000 rows the two cross near
  100 words at 16 members. A limit of 4 would send spans of 16 to 31 at 8000 rows to `std::sort`, 275
  against 215. The limit is interpolated, not timed at its edge; about one sort a sweep at n 8000 turns on it.

Candidates:

| candidate | effect | why bits cannot move, or why not taken |
|---|---|---|
| bitmap walk (taken) | the table; in place on the box, Q/P 0.941 warm and 0.934 cold | the one ascending arrangement |
| insertion sort as the fallback | beats libstdc++'s `std::sort` on the box from 16 members (202 against 275, 535 against 700), equals it below; loses to libc++'s on the Mac (572 against 534) | same array. Not taken: under the predicate a span of 16 or more falls back only above 384 rows a member, under one sort a sweep at n 8000, and a quadratic sort needs its own ceiling for caps above 64 |
| two-level walk (a summary bit a word) | Mac only (the critic's sortbench2): at 8000 rows 60, 111, 215 against the walk's 151, 230, 399; at 500 rows slower, 86, 138, 217 against 75, 107, 154 | same array. Not taken: no box timing, a second bitmap, and it wins where gp computations are rare (26 a sweep at n 2000, under 5 at n 8000: at most 0.4 percent of a sweep). A lever on a miss |
| cheaper in-order test | none: `std::is_sorted` leaves at the first descent, and 4 percent of spans are in order | same array |
| sort once a leaf a sweep (a memo by node, or the cache entry's list) | at most the draw's sorts, 39 percent | exact only with a membership check, itself a pass over the members against a table of n, or a hook in the moves |
| members kept in order at birth, change and swap (a stable partition) | removes the sort | NOT neutral: every leaf model sums in span order ([`Tree::partitionByPredicate`](../../src/bartcore/tree.hpp)); another tier, and the maintainer's question |

Prototype against the plain build on the Mac: nine gp fits `identical()` in draws, test fits, a copy's
continuation and the serialized state (scratch/gpsp/bits.R: 200 to 5000 rows; caps 64, 100 and 256; zero
weights; probit on two chains; three designated columns).

## Steps

1. model.hpp, ~50 lines: the predicate and its two constants, the two helpers with their asserts,
   `memberBits_` and its conditional sizing, the dispatch, the friend line; `sortedMembers`' comment says
   both methods return the one ascending arrangement.
2. [test_peer.hpp](../../tests/cpp/test_peer.hpp), ~30 lines: `TestPeer` reaches `sortedMembers`, each
   helper, `sortsByBitmap`, the size of `memberBits_`, and `holdsKernelEntry(leaf, tree, nodeIndex)` (the
   tree's entry at that node has members).
3. [test_model.cpp](../../tests/cpp/test_model.cpp), ~140 lines:
   - new `testGPLeafSortedMembers`. Leaves over 64, 110, 500, 4097 and 30000 rows at cap 64, and 500 rows
     at cap 256. EACH helper called directly at every span size from 2 to 64 (and 65, 128 and 200 at cap
     256; the bitmap helper wherever `memberBits_` is allocated, below the floor included), on spans in
     random order from a fixed generator that hold rows 0 and n - 1 and rows whose index mod 64 is 32 or
     more: every result equals a `std::sort` of the span. Two different spans in a row through the bitmap
     helper (a word left set would show). `sortedMembers` on the same spans, and on an in-order span, which
     returns the span's own pointer. The predicate's table at 15, 16 and 17 members and at 6 words a member
     and one more. `memberBits_` empty at 30000 rows and sized at the others. A leaf reinitialized from 100
     rows to 300 sorts a span of rows above 200.
   - [`testGPLeafMemberOrder`](../../tests/cpp/test_model.cpp), at 48 members: the score builds the entry,
     the span is reshuffled, `beginTreeDraw` runs, and `holdsKernelEntry` is read THEN, before any score or
     draw of the node (the next lookup rebuilds a dropped entry and would hide it). This pins without a
     clock the mutant the member-order review caught by timing alone, eviction against the unsorted span.
4. Records, ~25 lines: a dated paragraph in [gp-leaves.md](../design/gp-leaves.md) under
   [Stage 4 addendum: members in observation order (2026-10-09)](../design/gp-leaves.md#stage-4-addendum-members-in-observation-order-2026-10-09)
   (the two methods, the predicate's source, the box figures); TODO's item; this file's Landing note; the
   ledger entry. No NEWS and no help: gp leaves are new in 1.0-0 and nothing a user reads changes.

Reviewer's mutants, each to fail the named check: the word not zeroed after its walk (two spans in a row);
the shifted `1` left 32 bits wide (rows with index mod 64 of 32 or more; the undefined-behavior sanitizer
too); `& 63` or `* 64` wrong, or bits read from the top (every equality); the word count without its
round-up (row n - 1 at 110 and 500 rows, the helper's assert and the address sanitizer); `memberBits_`
not resized (the reinitialized leaf, same); the bitmap allocated where it cannot be used, or not where it
can (its size); a constant of the predicate moved (the table); the in-order return dropped (the pointer);
eviction against the unsorted span (the held entry). One mutant no test sees: `sortedMembers` calling the
comparison helper whatever the predicate says returns the same array, and only the box's reading below
catches it; the reviewer runs that reading on the mutant, or says it was not run.

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
  benchmarks/baselines/equivalence-308c45be.rds --bitwise --strict-coverage` 55 of 55, where
  ["GP leaves with a chi hyperprior"](../../benchmarks/R/equivalence.R) and
  ["GP leaves under NON-UNIT weights"](../../benchmarks/R/equivalence.R) are the two that build gp leaves
  (300 rows, cap 100: the bitmap from 16 members, `std::sort` below); `bcf-equivalence.R compare
  benchmarks/baselines/bcf-equivalence-1b7d730c.rds --bitwise` 15 of 15 and `multinomial-equivalence.R
  compare benchmarks/baselines/multinomial-equivalence-80b1c8d4.rds --bitwise` 11 of 11 (no gp leaf); the
  four `test-reproducibility-*.R` files unchanged (no gp leaf). Nothing is re-recorded.
- scratch/gpsp/bits.R on the base and on the slice, with a tenth fit at 30000 rows and cap 64 (no bitmap):
  ten of ten `identical()`.
- `R CMD check --as-cran` from a clean tarball; lintr, air, rc-codoc, win-drift, doc-freshness.

Timing, on the quiet box, orchestrator-run, shipped builds each `--preclean` into its own library; one
thread, one chain; harness scratch/gpsp-ab4.R (ab2.R with a seed range and the reps as arguments), arm
`cap64-n500-t50`, seeds `25:48`, 3 reps:

1. P' (the commit the slice branches from) against the slice tip. The states are burned in on P', so the
   chains are identical (the csv's two share columns agree in every row) and the ratio is the code alone.
   Bound: pooled warm at or under 0.98 (the prototype read 0.941, se 0.003 across seeds). The only bound.
2. A = da19dec4 against the tip, reported beside it with its per-seed range and no bound (B' read 1.017 on
   these seeds; about 0.96 is expected).

If 1 misses, run sortbench on the box at the tip's constants and report before changing anything.

## Budget and stops

Planned ~245 lines; recent engine slices ran about twice their plans, mostly in tests, so forecast 370 to
490. Implementer about three hours; gates about three hours of machine time; the box about twenty minutes.

Stop and report when: the diff passes ~500 lines; an equivalence scenario or a snapshot moves, or bits.R
differs (the invariant is broken, not a baseline stale); a path is found that hands a gp leaf a repeated
member or one at or above `numObservations_`; a span over the cap reaches `sortedMembers`; the change needs
anything outside `GPGaussianLeaf` and the tests; reading 1 misses.

## Calls made in planning

- **A floor of 16 members and 6 words a member, one rule.** The alternatives: the prototype's 8 and 4,
  which loses at spans of 8 to 15 from 2000 rows; a second tier taking spans of 8 to 15 where the bitmap
  has at most a word a member, which keeps the 20 ns at 500 rows. The cost is that 20 ns a span and a limit
  interpolated between measured rows.
- **`std::sort` stays the fallback, and the walk has one level** (the candidates table). The cost is up to
  a quarter of the fallback's time on libstdc++ where it runs, and whatever a second level would save from
  2000 rows up; neither was timed in place.
- **The bound is the identical-chains reading on seeds 25 to 48, alone.** The alternative keeps a bound on
  the tip against A, whose per-seed range on these seeds is 0.916 to 1.172 for B', or bounds on the other
  six arms, which fail by noise. The cost is that the slice is not held to the 1.05 it was opened on; that
  reading is reported, and B' already sits under it over 48 seeds.
- **The bitmap is not allocated where the predicate cannot pass.** The cost is a condition in
  `reinitialize` and an assert in the helper; otherwise one bit an observation per gp forest per chain.
- **Test hooks through `TestPeer`, not a counter in the leaf.** The alternative carries a kernel-build
  count beside the fallback tally, state no production path reads. The cost is a friend line.
- **bench-sampler.R compare is waived**, as it was for the member-order slice (dec-A192): none of its arms
  builds a gp leaf. The cost is that the hot-path rule's own gate does not run; reading 1 does.

## Not verified

- The predicate as chosen has not run anywhere: the box timed the prototype's rule in place and the two
  methods apart. At 500 rows the two rules differ only on spans of 8 to 15.
- The protocol's swapped-order run of Q against P; the two-level walk on the box; insertion sort in place.
- The prototype ran bits.R and the harness only: not tests/cpp, the tinytest suite, the equivalence trio,
  the snapshots or a sanitizer. The helpers, the asserts and the conditional allocation are unwritten.
- That no path hands a gp leaf a repeated or out-of-range member: read from the partition code and the
  leaf's own indexing of `u_`, not tested; the helper's assert is the first check of it.
- The Mac's times (the 2.5 percent, the 6 us, the two-level figures) were taken under load. The revision
  has had no second critique.
