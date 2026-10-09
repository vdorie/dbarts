# gp-copy-continuation: gp leaves compute over members in observation order

Status: PLANNED 2026-10-09.

agent: opus implementer, one; blind critique of this plan first; one opus reviewer who runs the mutants below.
rng: SHIFTING for every fit with gp leaves: a leaf under the size cap whose members the tree holds out of
observation order (in practice every gp fit after its first accepted birth) draws different values from the
same distribution. NEUTRAL, bit for bit, for constant, linear, vector and multinomial leaves (the change is
inside `GPGaussianLeaf` only) and for gp leaves whose spans are already sorted.
window: before the merge; engine slice after probit-k-mixing; engine slices stay serial.
budget: ~400 lines (model.hpp ~110, tests/cpp ~140, tinytest ~90, comments and docs ~50, MANIFEST ~10).

## Goal

A copy, a setState and a reload of a sampler with gp leaves continue the source to the last few digits, as
the copy and state help already promise. Every gp computation on a leaf under the cap reads the leaf's
members sorted by observation index, so a leaf's score, draws and test fits depend on its membership and not
on the history-dependent order of the tree's span. No state format, C API or facade change. The tier is
"Changes draws" ([Process by risk](README.md#process-by-risk)).

## Context

- Cause (diagnosis in scratch/gpcopy/diagnosis.md): a restore rebuilds each live tree from the identity order
  ([`Chain::rebuildLiveForest`](../../src/bartcore/chain.hpp), [`Tree::initialize`](../../src/bartcore/tree.hpp),
  [`Tree::repartitionSubtree`](../../src/bartcore/tree.hpp)), while the live order comes from
  [`Tree::partitionByPredicate`](../../src/bartcore/tree.hpp)'s unstable swap partition.
  [`GPGaussianLeaf::drawFromPosteriorForNode`](../../src/bartcore/model.hpp) builds the kernel and its
  Cholesky factor in member order and assigns standard normals by position, so a permuted span is a different
  valid draw. Everything else a restore carries (generators, lengthscales, standardization, k, sigma) is exact:
  handing the copy the source's span order alone makes it continue to 2.4e-15.
- Measured on bartcore e9868ce9 (scratch/gpcopy/paths.R, paths2.R): copy, setState and reload differ from an
  uninterrupted twin by 0.61 over 10 draws; with test rows, zero weights, probit and two chains, 0.56 to 3.4.
  Constant leaves differ by 1.8e-15.
- Prototype (scratch/gpcopy/protopkg, `// PROTO` lines in its model.hpp) sorts the tree's span in place for
  the length of each call and restores it after. It continues to 6.4e-15 (paths.R) and, in paths2.R, 6e-15 to
  1.9e-13 on train and 5e-12 to 1.6e-11 on test rows (alpha = C^-1 f amplifies rounding through the nugget's
  conditioning). test-gp-leaves.R and tests/cpp pass on it.
- Prototype speed, re-measured 2026-10-09 (loaded machine, interleaved, two rounds, min of 3; speed.R, n = 2000,
  50 trees, 77 to 92 percent of evaluations over the cap): base 0.074 / 4.10 s at caps 64 / 256, prototype
  0.085 / 5.01 s. Most of the loss is the kernel cache, not the sort:
  [`GPGaussianLeaf::evictStaleKernelEntries`](../../src/bartcore/model.hpp) runs in `beginTreeDraw` outside the
  prototype's sort and compares each entry's (sorted) members with the tree's span, so it evicts every cached
  leaf not already in order and the draw rebuilds its kernel and factor each sweep. With that comparison made
  in sorted order the prototype runs 0.080 / 4.06 to 4.20 s; the cap-64 remainder is presumably its
  per-call heap allocation and span restore (not measured apart).
- Alternatives not taken: storing member order in the state (option A in the diagnosis; a format change and
  about 50 percent more gp state), and sorting the tree's span permanently (a leaf model would mutate a const
  tree, and over-cap and constant-fallback sums would change order too).

## Constraints

- No state format, C API ([dbarts.h](../../inst/include/dbarts/dbarts.h)), facade virtual or R surface change.
  The public entry points of `GPGaussianLeaf` and the
  [`FunctionLeafModel`](../../src/bartcore/model.hpp) concept keep their signatures.
- States and saved trees written before the change install and replay unchanged: fits are stored per
  observation, and each saved leaf block carries its own alpha and rows in a matching order.
- Over the cap nothing changes: the constant fallback reads the tree's span as today and pays no sort.
- No allocation per call in steady state: the scratch keeps its capacity, as `cholV_` and the rest do.
- Out of scope: [`Chain::functionLeafValues`](../../src/bartcore/chain.hpp) sums a leaf's fits in span order
  for the flat record's per-leaf mean; a copy differs there by rounding only, inside the tolerance below.
- No NEWS: gp leaves are new in 1.0-0. No help change: the copy, setState and state items already say "to
  the last few digits, not bitwise", which this makes true for gp.

## Change

All in [`GPGaussianLeaf`](../../src/bartcore/model.hpp).

1. Sorted scratch. A private `mutable std::vector<index_t> memberScratch_` and a helper
   `const index_t* sortedMembers(const Tree&, const Node&, std::size_t numObs) const` that copies the span into
   it and `std::sort`s it, returning its data. It grows on demand and keeps its capacity (reserved to
   min(maxLeafSize_, n) in `initialize`), so after the first sweep it never allocates. One live pointer at a
   time: no caller nests two calls.
2. Entry points, each after its empty and over-cap branches and before anything reads members:
   - [`GPGaussianLeaf::logIntegratedLikelihoodForNode`](../../src/bartcore/model.hpp) (score);
   - [`GPGaussianLeaf::drawFromPosteriorForNode`](../../src/bartcore/model.hpp) (posterior draw);
   - [`GPGaussianLeaf::drawFromPriorForNode`](../../src/bartcore/model.hpp) (prior draw);
   - [`GPGaussianLeaf::appendLeafBlock`](../../src/bartcore/model.hpp) (flatten and live predict between runs).

   Each passes the sorted pointer down. `fits[members[r]]` writes are unchanged in meaning.
3. Helpers take `const index_t* members` in place of `(tree, node)`:
   [`GPGaussianLeaf::gatherLeafCovariates`](../../src/bartcore/model.hpp),
   [`GPGaussianLeaf::appendLeafRows`](../../src/bartcore/model.hpp),
   [`GPGaussianLeaf::cachedKernelForNode`](../../src/bartcore/model.hpp) (the entry key is the sorted list),
   [`GPGaussianLeaf::kernelAndFactorForNode`](../../src/bartcore/model.hpp),
   [`GPGaussianLeaf::logIntegratedLikelihoodOverPositiveWeights`](../../src/bartcore/model.hpp) and
   [`GPGaussianLeaf::drawFromPosteriorOverPositiveWeights`](../../src/bartcore/model.hpp) (`positiveScratch_`
   offsets index the sorted list). `anyZeroWeight` is order-free and may keep reading the span. Keep the two
   lines mutation-battery.R's ["poison 16"](../../benchmarks/R/mutation-battery.R) matches
   (`double w = ... weights[i];` then `double noise = residualVariance / w;`) textually intact, or re-match m16.
4. Draw cache. A `mutable std::vector<index_t> memberBuffer_` parallel to `alphaBuffer_`:
   [`GPGaussianLeaf::cacheAlphaForNode`](../../src/bartcore/model.hpp) takes the sorted pointer and appends it at
   the alpha's offset; [`GPGaussianLeaf::beginTreeDraw`](../../src/bartcore/model.hpp) clears both, capacity
   kept. [`GPGaussianLeaf::fitForTestObservationForNode`](../../src/bartcore/model.hpp) and
   [`GPGaussianLeaf::appendLeafBlockFromCache`](../../src/bartcore/model.hpp) read members from `memberBuffer_`,
   never from the tree. The test-row pool reads both buffers concurrently after the tree's draws, as it reads
   `alphaBuffer_` today; nothing writes them meanwhile.
5. Eviction. `evictStaleKernelEntries` compares an entry with the node's members in sorted order (sort the span
   into `memberScratch_`; entries exist only for leaves of 32 to cap members). This is the prototype's speed
   defect and the one place a missed sort costs time without changing a value.
6. Comments: the class comment gains one paragraph on member order (computations under the cap read members
   by observation index, so a leaf's draw is a function of its membership and a restore continues it);
   `drawFromPosteriorForNode` ("row order"), `appendLeafBlockFromCache` ("member order"),
   `fitForTestObservationForNode` ("the member ordering at draw time"), and `CachedLeafKernel` ("membership
   order") are restated.

Draw cache invariant: for a node drawn under the cap, `nodeAlphaOffset_[i] = o >= 0`, and for r in [0, m)
`alphaBuffer_[o + r]` is the weight of observation `memberBuffer_[o + r]`, the m members in ascending index,
with alpha = C^-1 f for the C built over that order. A constant node has offset -1 and no members. The saved
block's rows ([`addFlatFunctionPredictionsBelow`](../../src/bartcore/tree.hpp) replays them) follow the same
order, so replays still bit-match the live test fits.

## Tests

tests/cpp ([test_model.cpp](../../tests/cpp/test_model.cpp)):
- New `testGPLeafMemberOrder`: two leaves initialized alike over trees holding the same members, one in
  identity order and one shuffled, with test rows built as in
  [`testGPLeafFormats`](../../tests/cpp/test_model.cpp). At 8 members (no kernel cache) and 48 (cache
  engaged, cap 60), unit weights, positive weights and weights with zeros: score bitwise equal; posterior and
  prior draws from identically seeded generators give bitwise equal fits per observation and equal draw stats;
  every test row's fit bitwise equal; `appendLeafBlockFromCache` and `appendLeafBlock` blocks bitwise equal. A
  second score on the shuffled tree after reshuffling the span serves the cached kernel bitwise. Over the cap:
  one normal consumed, fits uniform.
- [`testGPLeafKernelCache`](../../tests/cpp/test_model.cpp): its state round trip now asserts continuation.
  After `cold.setState(state)` both samplers run 5 sweeps; fits and sigma agree to 1e-12. Its header comment and
  the "Only the install is asserted" comment are rewritten: a restored clone now draws over the same order.

tinytest:
- [test-gp-leaves.R](../../inst/tinytest/test-gp-leaves.R), a new section from paths.R and paths2.R, in
  process: n = 150, 15 trees, test rows, `n.chains = 2L, n.threads = 2L`; cases plain, weights with zeros, and
  probit. For each, an uninterrupted twin runs 3 + 10 draws; copy (with `updateState` TRUE and FALSE), setState
  into a sampler of another seed, and saveRDS/readRDS after `storeState` each continue it: train and test
  within 1e-8 absolute over the 10 draws (measured at most 1.6e-11; before the fix at least 0.56). The zero-weight
  case's two known warnings are expected by pattern and counted
  ([Gate hygiene](README.md#gate-hygiene)), not suppressed. Budget about 2 s.
- [test-mutate-then-serialize.R](../../inst/tinytest/test-mutate-then-serialize.R): the continuation check
  ["given the carried rng the continuation is the live one"](../../inst/tinytest/test-mutate-then-serialize.R)
  runs for gp as well as linear (the `if` goes).

## Baselines

- Moves: [equivalence.R](../../benchmarks/R/equivalence.R)'s two gp scenarios,
  ["GP leaves with a chi hyperprior"](../../benchmarks/R/equivalence.R) (gp) and
  ["GP leaves under NON-UNIT weights"](../../benchmarks/R/equivalence.R) (wtgp), both in the current
  equivalence-deb3fe50 ([MANIFEST](../../benchmarks/baselines/MANIFEST)). The diagnosis's "no snapshot pins gp"
  holds for the snapshots, not for this baseline.
- Does not move: the other 53 equivalence scenarios, bcf-equivalence-1b7d730c (15) and
  multinomial-equivalence-80b1c8d4 (11), none of which builds a gp leaf; the four test-reproducibility files
  (no gp); every exact gate (none builds a gp leaf). Any other mover is a defect: stop.
- Re-record: gp and wtgp fresh on the reference build (`EQUIVALENCE_SCENARIOS=gp,wtgp`), merged into a copy of
  deb3fe50 in its scenario order, named after the slice's code commit; deb3fe50 demoted to historical. Partition
  against deb3fe50 in z mode: 53 of 55 identical, gp and wtgp with no |z| above 4. The new file reproduces 55
  of 55 under `--bitwise --strict-coverage` from a second `--preclean` install.
- Oracle (MANIFEST rule P17): sbc.R's gp, gp-weighted and gp-mixed arms at the
  [Tier C - GP leaf](sbc-calibration.md#tier-c---gp-leaf) settings, each passing its ecdf band (about 23, 23 and
  6 min single-threaded; two at a time). Beside it the identity the MANIFEST row states: the sort is a
  permutation fixed by membership alone, and the score and both draws are equivariant under a fixed
  permutation, so the stationary distribution is unchanged.
- mutation-battery.R m16 still KILLs against the new baseline.

## Speed A/B

Same machine, quiet (maintainer-run or a quiet grant), shipped builds of the base (the slice's merge base) and
the slice tip, each `--preclean` into its own library, one thread, one chain.
- Grid: cap 64 and 256; n 500, 2000 and 8000; 10 and 50 trees; one designated column, plus one arm at three
  designated columns (n 2000, cap 256, 10 trees). 50 warm-up sweeps (caches warm), then 200 timed.
- Draws differ between the builds, so each arm runs three seeds and the statistic is the median over seeds and
  five interleaved repetitions (ABAB, both orders) of slice / base.
- Each arm reports its gp share (1 minus the fallback share, `dbartsGPFallbackWarning`'s tally); each cap needs
  at least one arm with a gp share of 0.5 or more, or the grid gains one.
- Accept: every arm at 1.03 or less. Levers before stopping: `std::is_sorted` on the span to skip the copy and
  sort; at eviction, keep an entry on size and bottom-ness alone (a lookup re-validates it anyway).
- bench-sampler.R compare is owed by the contract for a hot-path change, but none of its arms builds a gp
  leaf; this A/B replaces it (Open calls).

## Docs

- [gp-leaves.md](../design/gp-leaves.md): a dated note after
  [Stage 4 addendum: the kernel-cache pin's vehicle (2026-09-10)](../design/gp-leaves.md#stage-4-addendum-the-kernel-cache-pins-vehicle-2026-09-10)
  says gp computations under the cap read members by observation index, so a restored clone continues and its
  kernel matches; the addendum's "depends on the member list's ORDER" and closing paragraph point to it.
  [Stage 4 landing notes: kernel caching (2026-07-05)](../design/gp-leaves.md#stage-4-landing-notes-kernel-caching-2026-07-05)
  ("values and order"), [Stage 1 landing notes (2026-07-04)](../design/gp-leaves.md#stage-1-landing-notes-2026-07-04)
  ("f0's in row order first") and [Stage 2 landing notes (2026-07-04)](../design/gp-leaves.md#stage-2-landing-notes-2026-07-04)
  ("member order") are restated.
- The `testGPLeafKernelCache` comments (Tests). TODO's gp-copy-continuation item goes at landing.
- None of the help changes (Constraints).

## Steps

1. Change 1 to 6, with `testGPLeafMemberOrder` and the `testGPLeafKernelCache` assertion; tests/cpp green;
   `R CMD INSTALL --preclean -l <lib> .`; the full tinytest suite green.
2. The tinytest section and the test-mutate-then-serialize.R change; mutants below each fail a test.
3. Comments and gp-leaves.md.
4. After review: the re-record, MANIFEST row and oracle runs, in their own commit.

## Gates

On the slice tip against its own library, independently of the implementer
([RNG classes and their gates](README.md#rng-classes-and-their-gates), shifting):
- tests/cpp plain and under `-fsanitize=address,undefined`; the R-loaded ASAN path over test-gp-leaves.R and
  test-mutate-then-serialize.R (new buffer indexing).
- Full tinytest suite; `R CMD check --as-cran` from a clean tarball; lintr, air, rc-codoc, win-drift,
  doc-freshness.
- Reference build: the equivalence trio as in Baselines; the four snapshot files unchanged.
- The SBC oracle and the speed A/B above.

Reviewer's mutants, each of which must fail:
- the sort dropped from each entry point in turn (score, posterior draw, prior draw, appendLeafBlock);
- `fitForTestObservationForNode` or `appendLeafBlockFromCache` reading the tree's span instead of
  `memberBuffer_`;
- the positive-weight paths indexing the unsorted span;
- eviction comparing against the unsorted span: values unchanged, so only the cap-256 A/B arm catches it;
  the reviewer runs that arm on the mutant to show the A/B discriminates.

## Stop conditions

Stop and report when: the diff passes ~700 lines; a scenario other than gp and wtgp, or a snapshot, moves;
any A/B arm stays above 1.05 after the levers; an SBC arm flags on two seeds; the fix needs a state, facade or
C API change.

## Interactions

- probit-k-mixing (before this): changes probit draws and re-records the equivalence baseline. Disjoint code
  (this slice touches only `GPGaussianLeaf` and tests); whichever lands second re-records against the other's
  file, and this slice's partition is then against that file.
- state-install-keeps-spread: install paths only; no overlap.
- suite-stray-warnings: the new section adds no unpinned warning.

## Open calls

- Waive bench-sampler.R compare for this slice in favour of the gp A/B (no bench arm reaches the change)?
- SBC gp arms as the P17 oracle (about 50 min on two cores), or name the permutation identity alone?

## Estimate

Implementer about half a day; gates about four hours of machine time, an hour of it on a quiet machine.
