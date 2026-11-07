# Review 4, slice S-A: engine robustness

Status: built, in review (branch review4-engine-robustness).

rng: neutral, with one exception. Every complete run, and every state the package writes, draws bit for
bit as before; the equivalence baselines and the four snapshot files are unchanged ([Gates](#gates)).
The exception is one warm-start route: a donor whose DART weights hold 0 on a column the recipient's
forest may split on, as a forest under a column restriction writes them, now starts the recipient's DART
weights fresh (1/3 each over three columns, where main keeps the donor's 0.5, 0.5, 0), so the
draws after such a warm start differ ([Calls made](#calls-made), 3). What else changes is what an
incomplete run leaves, what an allocation failure leaves, and which edited or damaged states and
warm-start donors are refused.

Findings fixed: engine-01, engine-02, engine-03, engine-V-01, bridge-01 and bridge-V-01 of the fourth
whole-repo review, as amended by its blind critique (amendments 1 to 12 binding; 13 to 16 taken or left in
[Calls made](#calls-made)). Rulings this rests on: dec-B283 (a restore refuses a state whose values are
plainly broken), dec-B581 (a stopped run's results and saved trees are discarded), dec-A174, dec-A220 and
dec-A223.

## Fixes

### engine-01: a stopped run drops only the saved trees it wrote

Before: a run that stopped early (an interrupt, a draw callback's nonzero return, a sweep closure's stop,
an exception) returned without moving the sample cursor or the recorded-draw count, while the saved-tree
ring already held the slots it had rewritten, so `predict` and `getTrees` served those as earlier draws.

Now: `Chain::keptDrawsThisRun` counts a ring slot at the top of every kept sweep, after the cancel and
hook checks and before the sweep's first draw, so a sweep stopped part way (a monotone count, a
multinomial glue draw, an allocation failure between saved-tree writes) counts the slot it began. A guard
in `Sampler::run` applies `recordedDraws = K >= C ? 0 : min(R, C - K)` on every exit but the completed one,
with K the most any chain began, C the capacity and R the count before the run, the cursor unchanged. The
store then serves the earlier runs' draws in the slots no chain of the stopped run touched, oldest first.
The same guard clears a running flag: `Sampler::run` throws `std::logic_error` when it is called from
inside a run's sweep hook (amendment 13).

Tests: `testStoppedRunKeepsEarlierDraws` (test_sampler.cpp) - two chains on two threads with chain 0
stopping at draw 3, asserting the bound and that every served draw is one of run 1's in order (amendment
11c); a burn-in stop leaves the store unchanged; a stop at kept sweep 5 keeps C - 5 draws; a half-full
ring and a stop at draw 2 leaves it unchanged, and a full run afterwards writes from the unmoved cursor;
a run started in the hook throws and the outer run completes. `testMonotoneCountFailureKeptDraws`
(test_monotone.cpp) - an allocation failure inside a monotone count part way through a kept sweep, a full
ring of 6: 5 draws served, each run 1's next (amendment 11b). test-interrupt-saved-trees.R - the sweep
closure stopping before kept sweep 9 serves run 1's draws 10 to 20 through `predict` and `getTrees`, a
burn-in stop changes nothing, `storeState`, `saveRDS` and `readRDS` after the stop give the live
sampler's predictions (amendment 11a), a half ring is kept. test-capi.R - the compiled callback stopping
a full store of 10 at its draw 5 leaves the 4 draws it did not reach.

### engine-02: the sweep-closure entry writes only into buffers that fit

`bartcore_runWithCallback` checks each non-null result buffer before the run (`checkSweepBuffer`): sigma
and k one double a draw, train and test n or nTest times the reported locations, varcount numPredictors
integers, varprobs numPredictors doubles; a test buffer given with no test rows is refused the same way,
`results$<name> must be a <type> vector of length ...`. The closure's hook compares the sampler's shape
after each call; a change of n or nTest stops the run before that sweep stores anything and the entry
raises `the sweep callback changed the number of observations or test observations, which a run's result
buffers are sized by; change them between runs`, tested before the interrupt re-raise (amendment 7).

Tests: test-sweep-order.R - a wrong length or type for each of the six channels refused by name and the
sampler then running as an untouched one; correct buffers filled; a closure `setData` to another n and a
`setTestPredictor` to another nTest each raise the message, and the sampler runs again at the new shape;
a same-n `setData` inside the closure runs.

### engine-03 and engine-V-01: an allocation failure leaves totals equal to the trees

Before: `Chain::run`'s `bad_alloc` catch re-summed the totals for the monotone leaf only, so every other
leaf could leave a sweep with `totalFits` off its trees; a re-sum alone would read past a leaf table a
birth had outgrown (amendment 4).

Now: `Chain::reserveForMove` reserves, before each tree's move, room for one more child pair in the
tree's nodes, its leaf table and its bottom scratch, so nothing from an accepted move to the end of that
tree's draw reallocates them. `Chain::settleAfterAllocationFailure` (noexcept) makes the tree in flight
(`treeInFlight_`) consistent without allocating - a constant leaf's map, the monotone one's included,
marked stale and its table resized within the reserve, a vector leaf's block likewise - and re-sums
every forest through `rebuildTotalFitsFromTrees`, which gains the dense branch for linear and gp leaves.
`sampleNodeParametersFromPrior` catches for every leaf (amendment 10). `Tree::compactMaskPoolIfNeeded`
reserves its scratch (`liveMaskWordsBelow`) before pointing a rule into it.

The harness found two older defects of the same kind, fixed here: the gp kernel cache set an entry's
members (its key) before building its kernel, so a failure left an entry a later leaf matched with no
kernel (a segfault); `MonotoneLeafGeometry::build` left its axis map stale on a throw, which the next
build read past its axes (a segfault). The cache now writes the key last, and the geometry resets the map
on every exit. The fix round's harness found a third: the settle left the monotone leaf's map unmarked,
on the view that its failures put the tree back, but a failure in its branch score (the leaf list it
gathers) leaves the proposed tree in place, so the next sweep scattered sufficient statistics through a
map naming a split as a leaf (wrong statistics under NDEBUG, a failed assertion in tests/cpp). The
monotone tree's map is now marked stale too, and its next draw restarts its leaf values at the all-equal
zero, as a stale monotone tree does.

Tests: test_alloc.cpp (new suite `alloc`) replaces every global `operator new` and `delete` form, counts
the allocations of a sweep, and for each N fails the N-th on a fresh sampler, then checks the totals
equal the trees' sum, the state it leaves reads and installs on a new sampler, the failed sampler's next
sweep is consistent, and after it each tree's fits equal those of a new sampler built from its state, so
every leaf map the failure left behind its tree is rebuilt before it is read. The install goes to a new
sampler because an install on the failed one rebuilds its maps and hides a stale one. Arms: constant,
linear, gp and monotone leaves, each over a sweep of births and deaths from stumps and over
`sampleNodeParametersFromPrior`; a pooled-factor arm, a column of 60000 categories whose masks live in
the tree's pool, so a split on it crosses the pool's high-water mark and the failures land in the
compaction's reserve and around the stale-map marking; a tree-level arm compacting a pool of two live
masks stored past garbage, after each failure every rule reading its own words; and a kept-draws arm
whose failures land inside the kept sweep with a full ring. The sweep arms' birth probability is 0.5: at
1 a birth has no reverse death and is always rejected, so the trees never grew. All N fail cleanly:
constant 47 of 47 and 2 of 2, linear 85 of 85 (the prior draw makes no allocation), gp 102 of 102 and 58
of 58, monotone 189 of 189 and 105 of 105, pooled factor 75 of 75, compaction 1 of 1, kept draws 40 of
40 (39 inside the kept sweep).

### bridge-01 and bridge-V-01: values no chain can run from are refused

`Chain::drawnScalarsAreValid`, judged only where the chain draws the value and the state holds one: sigma
positive and finite; each forest's k positive and finite, except 0 where its prior lets a draw reach 0
(`Chain::kDrawCanBeZero`, a squared scale below 2^-960: the draw's rate overflows below 2^-1025, and
below 2^-969 a gamma draw can round to 0) and +Inf where it lets a draw reach +Inf
(`Chain::kDrawCanBeInfinite`, amendment 1); DART weights finite, not negative, not all zero, and positive
on every column the forest may split on, the masked columns free (amendment 3), the concentration
positive and finite; the response transform finite. `Chain::generatorStateIsValid` wraps
`ext_rng_serializedStateIsValid`: the Mersenne Twister index in [0, N] and not every essential bit zero
(amendment 14). Ordinal cutpoints: finite, the first exactly 0, strictly increasing (amendment 5).
`Chain::stateIsValid` applies all of them, `Sampler::setState` judges the generators before the
missing-direction coins are drawn from them, and setState, a copy and a reload refuse with `state is not
consistent with this sampler`. `ext_rng_readSerializedState` re-derives the gamma memo from the shapes it
records (`ext_rng_rederiveGammaMemo`), so an edited memo constant cannot reach a draw (amendment 6).

Warm start: `Sampler::installForests` holds each chain's assembled state to the same floors and returns
`WarmStartResult::valueMismatch`, which the bridge raises as `warm-start donor holds a value this sampler
cannot run from: ...; nothing was installed`, with and without force (amendment 12). A DART donor whose
well-formed weights hold 0 on a column this forest may split on, as a forest under a column restriction
writes them, starts the recipient's DART fresh (`Chain::startDartFresh`) instead.

Tests: `testStateValueFloors` (test_state.cpp) - each floor refused at its edge and its package-written
value accepted, in both forms; k = +Inf accepted under scales Inf and 1e150, refused under 1e140 and 2;
k = 0 written under 1e-160 and accepted by setState, a new sampler and a warm start, refused under 1e-140;
the index at -1, -100000, N + 1 and N + 2 refused, 0 and N accepted; all-zero and low-bits-only words
refused; a memo edit installed recomputed to the bit and running bitwise equal; the DART, ordinal and
warm-start cases. test-state-value-floors.R - the same through setState in both forms, two copies of one
sampler compared by `identical()` (amendment 8), a reload refused rather than a crash, installTrees
refused by its message, k = Inf under chi(1.5, Inf) installed, the k = 0 a chi(1.5, 1e-160) sampler writes
taken by setState, a copy, a reload and installTrees, a restricted DART donor starting a
recipient's weights at 1/3 each.

## Gates

Run on arm64 macOS against the branch's own libraries (shipped and reference builds, each installed with
`--preclean`).

- tests/cpp: every suite passes, plain and under `-fsanitize=address,undefined` (445 checks reported
  ok in each after the fix round, 443 before it, no sanitizer report).
  The new suite `alloc` reports the counts under engine-03 above.
- Fix round, shipped build: test-state-value-floors.R, test-interrupt-saved-trees.R, test-capi.R,
  test-sweep-order.R and test-interrupt-signal.R, 78, 9, 324, 26 and 2 results, 0 failed. The full
  tinytest, R-loaded ASan and R CMD check lines below are from the build, before the fix round.
- Full tinytest, one process, shipped build: 41276 results in 303 files, 0 failed.
- Equivalence trio on the reference build, `--bitwise` (and `--strict-coverage` on equivalence.R): 55,
  19 and 11 scenarios compared, 0 skipped, every one identical draws; the same again after the fix round.
- The four snapshot files on the reference build: 14, 3, 7 and 3 assertions, 0 failed (unchanged); the
  same again after the fix round.
- R-loaded ASan (the bridge compiled with `-fsanitize=address`, run through `bin/exec/R` with the
  runtime inserted): the four new or touched files and nine neighbors (state, warm start, callbacks,
  ordinal, DART, gp, monotone, serialization), 1286 results, 0 failed, no report. test-interrupt-signal.R
  is not counted: its child process cannot load an instrumented library without the runtime; it passes
  at home on the shipped build.
- R CMD check --as-cran on the built tarball: Status 1 NOTE, CRAN incoming feasibility only (an
  external site's expired certificate, cited from two help pages this slice does not touch, and the Date
  field).
- `lintr::lint_package()` 0 lints; `air format --check .`, check-rc-codoc.R (53 methods), check-win-drift.R,
  check-doc-freshness.R and `mutation-battery.R verify-anchors` all OK; after the fix round lintr, air,
  codoc, win-drift and doc-freshness again OK.
- Mutants, each built into a copy of the tree and run against tests/cpp:

  | mutant | killed by |
  |---|---|
  | the guard dismissed on the cancel path | three stopped-run checks in `testStoppedRunKeepsEarlierDraws` |
  | the slot counted after the saved writes (amendment 2) | `testMonotoneCountFailureKeptDraws` and the alloc suite's kept-draws arm |
  | `if constexpr (NormalizedLeafModel<L>)` restored in `Chain::run`'s catch | the alloc suite's constant, linear and gp sweep arms |
  | each floor of `drawnScalarsAreValid` removed: sigma; k negative; k infinite; transform; DART well formed; DART positive where the forest splits; DART concentration | its own check in `testStateValueFloors`, the two DART weight floors each by the weights check |
  | k = 0 accepted under every scale; k = 0 refused under every scale (`kDrawCanBeZero` false) | the refusal under 1e-140; the acceptance of a chi(1.5, 1e-160) sampler's own state |
  | the settle's stale mark removed (the reviewer's M6) | the engine's map assertion in the alloc suite's constant sweep arm; with `-DNDEBUG`, the pooled-factor arm's routed-fits check |
  | the monotone leaf excluded from the stale mark again | the engine's map assertion in the alloc suite's monotone sweep arm |
  | the compaction's reserve removed | the alloc suite's mask-pool compaction arm |
  | each ordinal floor removed: first 0; finite; increasing | the cutpoints check |
  | the warm start's `valueMismatch` removed | the warm-start donor check |
  | each generator floor removed: the index range; the essential bits; the memo re-derivation | the index, essential-bits and memo checks |

  The floor checks reuse one sampler for a run of edits, so once a mutant lets an edit in, the later
  checks on that sampler fail with it; every mutant whose check follows another's on one sampler was
  built alone.

## Calls made

1. The gamma memo is re-derived from the shapes it records on every serialized read
   (`ext_rng_rederiveGammaMemo`, through the same out-of-line helpers the gamma sampler calls, so the
   arithmetic is the sampler's own), not zeroed as amendment 6 proposed. Zeroing keeps the draws but
   makes a state installed and stored again differ from the one given until the next gamma draw, so an
   `identical()` round trip of a sampler's own state fails; a check that the memo equals its
   recomputation would refuse a state saved on a build that contracts the arithmetic differently.
2. k = +Inf is accepted wherever a draw can overflow, `scale^2 > 2^970` (`Chain::kDrawCanBeInfinite`),
   not only under an infinite scale as amendment 1 states: under chi(df, 1e150) the draw's bound
   can overflow too, and the amendment's test would refuse a state such a sampler writes. Its mirror
   came in the review's fix round: under chi(1.5, 1e-160) the draw's rate overflows and the chain
   writes k = 0, so k = 0 is accepted wherever a draw can reach it, a squared scale below 2^-960
   (`Chain::kDrawCanBeZero`: below 2^-969 the smallest gamma draw the generator makes, 2^-106, times the
   gamma scale rounds to 0), at every door the floor reaches. Alternative: refuse such scales when the
   prior is built, a change to what a user may pass.
3. A warm-start donor whose DART weights are well formed but hold 0 on a column this forest may split on
   starts the recipient's DART fresh; setState refuses the same weights. Such a donor is package-written
   (a `forest(x1 + x2)` fit on the same three predictors, which the tinytest builds), and a warm start
   crosses model differences by design. Alternative: refuse it as `valueMismatch`.
4. After an allocation failure only the tree in flight is settled - its map marked stale, its leaf
   table resized inside the reserved room - rather than every constant-leaf tree's map (amendment 4's
   text). Every other tree finished its move and draw, so its map and table agree, and the in-flight
   marker is needed anyway to find the table to resize. Alternatives: mark every map stale (a full
   rebuild at each tree's next draw), or reset a tree with a short table to a stump.
5. `Chain::reserveForMove` grows geometrically, so a move where the room is there costs one compare.
   Alternative: reserve exactly two more nodes each move, which reallocates at every birth.
6. Amendment 13 taken: `Sampler::run` throws `std::logic_error` while a run is in progress. A run
   re-entered from the sweep hook resets every chain's slot base mid-run and breaks the ring's order;
   the flag is one bool beside the guard. Alternative: leave re-entry undefined, since nothing calls
   the entry.
7. Amendment 14 taken: index N accepted and N + 1 refused, and the essential bits (the top bit of the
   first word and every later word) not all zero. Alternative: the plan's "every word 0", which admits
   a state whose only set bits are the first word's low 31, a state the twister treats as zero.
8. Amendment 15 taken as a sentence in bart-as-a-component.md section 5; the entry is unexported.
   Alternative: refuse `storeState` while a run is in progress, an R-surface change for an internal
   entry.
9. Amendment 16 left: at any throw point the variance forest's combined variance equals the product
   of its factors within rounding, and the next sweep's closing recompute restores it exactly; a
   recompute in the catch would buy exactness only until that sweep.
10. The warm-start message departs from amendment 12's text: sigma must be "positive and finite"
    (sigma = Inf is refused), and the DART weights "finite, not negative and not all zero", since a
    well-formed donor with a 0 where this forest may split starts fresh (call 3) rather than being
    refused.
11. `Sampler::setState` judges the generators before `completeMissingDirections` draws its coins from
    them, besides `Chain::stateIsValid`'s judgment, since those draws come first. Alternative: judge only
    in `stateIsValid`, after a draw from a generator it would refuse.
12. The two model.hpp fixes (the gp kernel cache's key written last, the monotone geometry's axis map
    reset on every exit) are outside the plan's list. The harness found both as segfaults after an
    allocation failure, of engine-03's kind, and the gp and monotone arms cannot pass without them.
    Alternative: leave them to a later slice and drop those arms.
13. The harness builds a fresh sampler for each N rather than restoring one chain state: scratch
    capacity survives `setState`, so after the counting run most allocations no longer happen and the
    N-th would not be reached. The monotone kept-draws test lives in test_monotone.cpp, beside the
    `failNextCount` scaffolding it uses, not test_sampler.cpp; the harness is its own suite.
14. bart-as-a-component.md now agrees with architecture.md that `bartcore_runWithCallback` is kept as
    the entry host packages build on, closing the disagreement the plan left open. Alternative: retire
    the entry, which nobody proposed.

## Design-doc register

- docs/architecture.md, "Callbacks and errors": the sweep-closure entry runs every R call under
  `R_ToplevelExec`, refuses a wrong buffer and stops on a shape change; a stopped run drops the saved trees
  it wrote and keeps the earlier draws; an allocation failure re-sums the fits before it leaves the engine.
- docs/design/bart-as-a-component.md, section 5: `bartcore_runWithCallback` is the entry host packages
  build on (agreeing with architecture.md); `R_ToplevelExec`; the buffer and shape rules, a same-shape
  `setData` legal; a mutator inside the closure under `updateState = TRUE` stores a mid-run state.
- docs/design/per-draw-callbacks.md: the stop paragraph describes the guard (`Chain::keptDrawsThisRun`,
  `Sampler::run`) in place of the inconsistent-after-abort property; `R_ToplevelExec`.
- inst/include/dbarts/dbarts.h (comments only, no ABI change): `dbarts_draw_callback`,
  `dbarts_sampler_run` and `dbarts_sampler_numSavedSamples` say the run drops the saved trees it wrote, and
  that an allocation failure leaves the chains runnable.
- man/dbartsSampler-class.Rd: `run`'s value (amendment 9), `setState`'s refusals, `installTrees`' DART
  fresh start and value refusal.
