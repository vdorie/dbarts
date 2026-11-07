# Review 4, slice S-A: engine robustness

Status: LANDED on main 2026-11-07 (branch review4-engine-robustness).

rng: neutral, with one exception. Every complete run, and every state the package writes, draws bit for
bit as before; the equivalence baselines and the four snapshot files are unchanged ([Gates](#gates)).
The exception is one warm-start route: a donor whose DART weights hold 0 on a column the recipient's
forest may split on, as a forest under a column restriction writes them, now starts the recipient's DART
weights fresh (1/3 each over three columns, where main keeps the donor's 0.5, 0.5, 0), so the
draws after such a warm start differ ([Calls made](#calls-made), 3). What else changes is what an
incomplete run leaves, what an allocation failure leaves, which edited or damaged states and warm-start
donors are refused, and three cases main treats otherwise:

- setState refuses a state whose DART weights hold 0 on a column this sampler's forest may split on, a
  restricted-DART sampler's own state given to an unrestricted sampler among them; main installs it
  ([Calls made](#calls-made), 15).
- setState and a warm start refuse a k of 0 or +Inf under a prior whose draw cannot reach it, a
  chi(1.25, 1e-160) or chi(1.25, Inf) sampler's own state given to a chi(1.25, 2) sampler among them;
  main installs both, and its k = Inf chain stays at Inf on every later sweep (call 2).
- An edited gamma memo in a generator state is corrected on read rather than refused; main draws from
  it as edited (call 1).

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
tree's draw reallocates them. `Chain::settleAfterAllocationFailure` (noexcept) makes the chain the one
its state describes, without allocating. The mean tree in flight (`treeInFlight_`) keeps its structure,
its leaf table covering its nodes inside the reserve (a new leaf at 0); a tree being regrown or installed
(`treeInFlightRegrown_`), or one whose table cannot cover it, goes back to a single leaf at 0. Its map is
then rewritten from its partition (`writeLeafOfBelow`) and a vector leaf's fits from its blocks
(`writeParameterFitsBelow`), the monotone leaf's block restarting at the all-equal zero. A variance tree
being drawn or installed (`varianceTreeInFlight_`) goes back to a single leaf, its factor 1, the variance
surface the product of the factors again. Every forest is then re-summed through `rebuildTotalFitsFromTrees`, which
gains the dense branch for linear and gp leaves. The settle runs from `Chain::run`'s catch and from every
other path that resets and rebuilds a tree: `growForestFromRoot` (reached by `$growFromRoot` and
`bart(n.grow.sweeps =)`), `sampleTreesFromPrior`, `sampleVarianceForestFromPrior`, and the installs
behind setState and a warm start (`rebuildLiveForest`, `rebuildVarianceForest`), on a throw or a refusal
part way. The paths that move the store under the trees settle too, each as the audit below says; the
second fix round's records said they re-route inside their transactions and reset none, which was
wrong: the forced refresh and the data replacement rebuilt every tree with nothing to settle them.
`sampleNodeParametersFromPrior` catches for every leaf (amendment 10). `Tree::compactMaskPoolIfNeeded`
reserves its scratch (`liveMaskWordsBelow`) before pointing a rule into it.

The harness found two older defects of the same kind, fixed here: the gp kernel cache set an entry's
members (its key) before building its kernel, so a failure left an entry a later leaf matched with no
kernel (a segfault); `MonotoneLeafGeometry::build` left its axis map stale on a throw, which the next
build read past its axes (a segfault). The cache now writes the key last, and the geometry resets the map
on every exit. The first fix round's harness found a third: the monotone leaf's map was left unmarked,
but a failure in its branch score leaves the proposed tree in place. The second fix round found two
more. `growForestFromRoot` reset and regrew each tree outside `treeInFlight_`, so after a failure there
the next sweep failed the map assertion, or under `-DNDEBUG` corrupted the heap; the prior draws and the
installs had the same gap. And the stale-map settle of the build, which kept a tree's old fits until its
next draw, left the state read just after a failure describing a different chain (constant, linear,
monotone and pooled-factor arms), and a mutant skipping the monotone table's resize went unseen. Hence the
settle above.

Tests: test_alloc.cpp (new suite `alloc`) replaces every global `operator new` and `delete` form, counts
the allocations of a step, and for each N fails the N-th on a fresh sampler. After each failure, before
any further sweep, the totals equal the trees' sum, every leaf table covers its tree, every monotone tree
is in order, and a new sampler given the state it reads has its fits and variance surface; that state
installs, the failed sampler's next sweep is consistent, and after it each tree's fits equal those of a
new sampler built from its state. The install goes to a new sampler because an install on the failed one
rebuilds its maps and hides a stale one. Arms: constant, linear, gp and monotone leaves, each over a sweep
of births and deaths from stumps, `sampleNodeParametersFromPrior`, `sampleTreesFromPrior` and a state's
install, and, for the constant and monotone leaves, `growFromRoot`; a variance forest over a sweep, its
prior draw, an install and a grow; a monotone settle arm, a tree's leaves put out of order and settled
in flight; a pooled-factor arm, a column of 60000 categories whose masks live in the tree's pool, so a
split on it crosses the pool's high-water mark and the failures land in the compaction's reserve; a
tree-level arm compacting a pool of two live masks stored past garbage, after each failure every rule
reading its own words; and a kept-draws arm whose failures land inside the kept sweep with a full ring.
The sweep arms' birth probability is 0.5: at 1 a birth has no reverse death and is always rejected, so
the trees never grew. All N fail cleanly (sweep, prior leaf, prior tree, install, grow): constant 47,
2, 16, 32, 62; linear 85, 0 (the prior draw makes no allocation), 17, 32; gp 102, 58, 13, 23; monotone
189, 105, 15, 82, 69; variance forest 72 in a sweep, 15 in its prior draw, 41 in an install, 60 in a
grow; pooled factor 75, compaction 1, kept draws 40 (39 inside the kept sweep).

The third fix round adds checks after each failure: every leaf holds a row; the store's grid is whole,
each column holding as many points as its count says, and its codes are those of the values held or
the values moved to; a linear or gp leaf's covariates are the store's values under the constants the
leaf holds; the test fits, the variance surface and the leaf's test covariates cover the store's test
rows; every chain is checked, not only the first; and, in new per-move arms, every tree splits as the
unfailed sweeps have it at the start or the end of one, so a proposal left in place shows. Its new
arms, for each leaf: three sweeps of each move kind alone (births and deaths, changes, swaps, perturbs,
rule Gibbs steps) from trees grown 25 sweeps; a warm start; an unforced setPredictor and a forced one
on grids derived again. For each leaf and a variance forest: a forced tree refresh (setCutPoints); a
data replacement over a reversed column, with and without test rows (a variance forest's without); a
test data replacement; and a state install that repairs, into a sampler whose columns are squeezed
into the middle tenth of their grids, so the stored trees' outer leaves hold no row. For the linear and
gp leaves, a standardization rederive over those squeezed columns, so the constants move. For a
constant leaf, a variance forest and two chains: a forced and an unforced setPredictor, each on the
grids held and on grids derived again; a forced one deriving grids longer than the 7 and 5 cuts
setCutPoints gave; a per-row update; a state install and a warm start across grids; and, for the two
chains, the forced refresh and the data replacement. All 114 arms pass under
`-fsanitize=address,undefined`, every failure reaching the caller and none absorbed.

#### Audit: every path that resets, rebuilds, re-routes or regrows a tree

| path | a throw part way is settled by | arm |
|---|---|---|
| a sweep's moves (birth and death, change, swap, perturb, rule Gibbs), `Chain::run` | the move taking its proposal back before it rethrows (the birth's pair released, the death's node put back, the changed, swapped or perturbed subtree restored, the incumbent rule put back); then the settle keeps the tree in flight, rewrites its map and fits and re-sums the totals | the sweep arms (births; every move), the per-move arms, each leaf; the variance, pooled-factor and kept-draws arms |
| leaf draws, in a sweep and `sampleNodeParametersFromPrior` | the settle, the tree in flight kept | the prior leaf draw arms |
| prior tree draws, `sampleTreesFromPrior` and the variance forest's | the tree in flight and every tree after it, in its forest, later forests and the variance forest, back to a single leaf (`markRebuildAt`, `restUnreached_`) | the prior tree draw arms, the variance prior draw |
| the grow from the root | the same | the grow arms (constant, monotone, variance) |
| a state install, forced or not (setState; `installLiveForestTrees`, `installVarianceTrees`), the forced one repairing a tree (`rebuildLiveForest`, `rebuildVarianceForest`) | the same | the state install arms; the repairing install arms, each leaf and variance |
| a warm start (`installForests`) | the chains it touched restored from the states taken before it | the warm start arms |
| an install across grids (setState or a warm start from another grid, `installAcrossGrids`) | the store's grid put back; once the install has begun, every chain's trees back to a single leaf; the refresh after it as the next row | the across-grid install and warm start arms, one and two chains |
| the forced tree refresh (`forceRefreshTrees`, behind setCutPoints, an install across grids and a forced setPredictor) | the trees re-routed before the throw kept, the one in flight and every later one reset (marks per tree); every later chain reset (`refreshEveryChain`); a throw while setCutPoints installs the grid resets every chain, each column's grid going in whole (`setCutPointsForColumn`) | the forced refresh arms, each leaf, variance, two chains |
| a forced setPredictor or updatePredictor | before the refresh, every chain reset, a grid derived again going in whole per column (`deriveNumericCuts`); in the refresh, as the forced refresh | the forced setPredictor arms, on the grids held, on new ones and on longer ones |
| an unforced setPredictor or updatePredictor (`runPredictorTransaction`) | the store's snapshot taken whole before anything moves (`WholeMatrixUpdate`, `SubsetUpdate::applied`); a throw before the fits are rebuilt puts the store back and re-routes the trees; one while they are rebuilt settles every tree in place (`settleEveryTreeInPlace`) | the setPredictor arms, on the grids held and on new ones |
| the per-row predictor update (its session's finalize) | the cells committed, every tree re-routed over them and settled in place | the per-row update arms |
| a data replacement (setData; `applyNewData`) | the store's replacement all or nothing (`ColumnStore::setData`); a throw after it drops the test rows and resets every chain not yet given the rows; within a chain, as the forced refresh, the leaf covariates gathered again from the store's values, the variance forest after | the data replacement arms, with and without test rows, each leaf, variance, two chains |
| a test data replacement (`setTestData`) | the test rows dropped, the store's test count back to 0, so no chain's test fits fall short of it; the training trees are not touched | the test data replacement arms, each leaf and variance |
| a standardization rederive (`rederiveStandardization`, behind setPredictor and `rederiveAfterInstall`) | the covariates gathered again under the constants the leaf holds, every tree's fits rewritten from its blocks, the totals re-summed | the rederive arms, linear and gp, over moved values |

Not settled, each listed in TODO: a data replacement with another row count, where the store, each
chain's per-row storage, its response and its leaf covariates are sized in turn, so a throw between them
leaves them sized for different counts and the sampler must not be used again; and a standardization
rederive's saved draws, or a linear leaf's in a data replacement, where one the rewrite had not reached
keeps its old basis. Each is reached from R only through an allocation failure, which R reports as an
error from the call. A failed test data replacement, or a data replacement with test rows, leaves the
sampler with none: its next run records no test fits until test rows are set again.

The reviewer's empty leaf, left by a failure in the constant sweep arm on the second round's code, is
reachable from R: an allocation failure inside a change, swap or perturb proposal, between its refresh
and its decision, under `$run()` or a fit. R reports the error; the sampler then kept the proposed tree,
whose leaf might hold no row. A user who ran on would see no further error, but the chain would run on
a tree the moves never admit: the empty leaf held at zero, not drawn, until a later move removed it;
a save and load collapsed it, and setState refused the state it stored unless forced. Each move now takes its proposal back before it rethrows, so the tree is
the one before the proposal.

The reviewer's Y1 mutant, the rebuilt tree's reset removed, survived the second round's harness; its
records called it killed. The mutant the records ran removed the reset line whole, so also the reset of
a tree whose table cannot cover it, which the monotone install arm caught. Y1 skips only the rebuilt
tree's reset, and leaves empty leaves in the gp install arm, which nothing then checked. It now dies,
with assertions on and under `-DNDEBUG`, by the state checks (every leaf holding a row among them) in
the forced refresh, data replacement, gp install and across-grid install arms.

The round's harness found four more gaps. The predictor transaction's revalidation allocated its
working state outside its rollback, so a failure there left the store moved under trees routed by the
old one; it is now allocated inside. The predictor arms' shift of 0.004 had moved no row across a cut,
the column's cuts lying midway between values 1/80 apart, so those arms had re-routed nothing; the
shift is now 0.009, and the failure above showed only once it was. (The data replacement over the
shifted values had the same blind spot, a grid derived again absorbing the shift, so it now reverses a
column; the rederive arm, over values whose constants it derived again unchanged, moved nothing, so it
now runs over squeezed columns; and no arm derived a grid longer than the one it replaced, the only
case where the derivation allocates after its first write, so one now does.) The store's test build sets its count before it allocates, and each chain sizes its test fits
to that count, so a failure there left the next sweep writing past them (a segfault under
AddressSanitizer); a throw now drops the test rows. And a data replacement that failed before the leaf
was reinitialized left a linear or gp leaf's covariates gathered from the values the store had
replaced; the catch now gathers them again.

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

- tests/cpp: every suite passes, plain and under `-fsanitize=address,undefined` (548 checks reported
  ok in each after the third fix round, 460 after the second, 445 after the first, 443 before it, no
  sanitizer report).
  The new suite `alloc` reports the counts under engine-03 above.
- Each fix round, shipped build: test-state-value-floors.R, test-interrupt-saved-trees.R,
  test-capi.R, test-sweep-order.R and test-interrupt-signal.R, 78, 9, 324, 26 and 2 results, 0 failed.
  The R-loaded ASan and R CMD check lines below are from the build, before the fix rounds.
- Full tinytest, one process, shipped build: 41281 results in 303 files, 0 failed, after the second and
  third fix rounds (41276 at the build).
- Equivalence trio on the reference build, `--bitwise` (and `--strict-coverage` on equivalence.R): 55,
  19 and 11 scenarios compared, 0 skipped, every one identical draws; the same again after each fix
  round.
- The four snapshot files on the reference build: 14, 3, 7 and 3 assertions, 0 failed (unchanged); the
  same again after each fix round.
- R-loaded ASan (the bridge compiled with `-fsanitize=address`, run through `bin/exec/R` with the
  runtime inserted): the four new or touched files and nine neighbors (state, warm start, callbacks,
  ordinal, DART, gp, monotone, serialization), 1286 results, 0 failed, no report. test-interrupt-signal.R
  is not counted: its child process cannot load an instrumented library without the runtime; it passes
  at home on the shipped build.
- R CMD check --as-cran on the built tarball: Status 1 NOTE, CRAN incoming feasibility only (an
  external site's expired certificate, cited from two help pages this slice does not touch, and the Date
  field).
- `lintr::lint_package()` 0 lints; `air format --check .`, check-rc-codoc.R (53 methods), check-win-drift.R,
  check-doc-freshness.R and `mutation-battery.R verify-anchors` all OK; after each fix round lintr, air,
  codoc, win-drift and doc-freshness again OK.
- Mutants, each built into a copy of the tree and run against tests/cpp:

  | mutant | killed by |
  |---|---|
  | the guard dismissed on the cancel path | three stopped-run checks in `testStoppedRunKeepsEarlierDraws` |
  | the slot counted after the saved writes (amendment 2) | `testMonotoneCountFailureKeptDraws` and the alloc suite's kept-draws arm |
  | `if constexpr (NormalizedLeafModel<L>)` restored in `Chain::run`'s catch | the alloc suite's constant, linear and gp sweep arms |
  | each floor of `drawnScalarsAreValid` removed: sigma; k negative; k infinite; transform; DART well formed; DART positive where the forest splits; DART concentration | its own check in `testStateValueFloors`, the two DART weight floors each by the weights check |
  | k = 0 accepted under every scale; k = 0 refused under every scale (`kDrawCanBeZero` false) | the refusal under 1e-140; the acceptance of a chi(1.5, 1e-160) sampler's own state |
  | the in-flight marker removed from the grow from the root | the engine's map assertion in the constant grow arm; with `-DNDEBUG`, a segfault there |
  | the marker removed from the prior tree draw; from the install | with `-DNDEBUG`, a trap in the constant prior tree draw arm; a segfault in the constant install arm |
  | the variance marker removed from the prior draw; from the install | with `-DNDEBUG`, the variance prior draw arm's state check; the variance install arm's state and routed checks |
  | the settle's map rewrite removed (the reviewer's M6 on this settle) | with `-DNDEBUG`, the constant sweep arm's state check, then a trap |
  | the leaf table's resize removed (the reviewer's X2) | the constant sweep arm's state check, then a segfault, with and without `-DNDEBUG` |
  | the monotone leaf excluded from the map rewrite (the reviewer's unmarked mutant on this settle) | the engine's map assertion in the monotone prior tree draw arm; with `-DNDEBUG`, a trap there |
  | the monotone block's zero restart removed | the monotone settle arm, with and without `-DNDEBUG` |
  | the vector leaf's fits rewrite removed | with `-DNDEBUG`, the linear sweep and prior tree draw arms' state checks |
  | the settle's reset line removed whole, every reset with it | the engine's map assertion; with `-DNDEBUG`, the constant install arm's state check, then a trap |
  | the variance tree's reset removed | with `-DNDEBUG`, the variance prior draw and install arms' state checks |
  | the compaction's reserve removed | the alloc suite's mask-pool compaction arm, with `-DNDEBUG` |
  | each ordinal floor removed: first 0; finite; increasing | the cutpoints check |
  | the warm start's `valueMismatch` removed | the warm-start donor check |
  | each generator floor removed: the index range; the essential bits; the memo re-derivation | the index, essential-bits and memo checks |

  Third fix round, each against the alloc suite with assertions on and under `-DNDEBUG`, killed the
  same way in both unless the row says otherwise:

  | mutant | killed by |
  |---|---|
  | the reviewer's Y1, the rebuilt tree's reset skipped | the forced refresh, data replacement, gp install and across-grid install arms' state checks |
  | each move's undo removed: birth, death, change, swap, perturb, rule Gibbs | the per-move arms' structure check (linear leaf; the swap's in the gp arm); the rule Gibbs undo by the constant arm's state check |
  | the forced refresh's settle removed; the data replacement's | the engine's totals assertion; with `-DNDEBUG`, the forced refresh and data replacement arms' totals checks |
  | the forced refresh's per-tree marks removed; the data replacement's first mark | the forced refresh arms' state checks; the linear data replacement arm's |
  | the trees after the one in flight left as they were (`restUnreached_`) | the forced refresh arms' state checks |
  | the data replacement's covariate regather removed | the linear data replacement arm's state check |
  | the unforced transaction's rollback removed; the revalidation's rollback; the snapshot's flag; the revalidation state allocated outside its rollback | the setPredictor arms' state checks (new grids; held grids; held grids); a segfault in the constant setPredictor arm |
  | the forced transaction's reset removed | the forced setPredictor arms' state checks |
  | the per-row session's settle removed | the variance forest's per-row update arm's state check |
  | a later chain's reset removed: in `refreshEveryChain`; after a chain's data replacement; in the install across grids; in setCutPoints; in setData before the chains | the two-chain growing-grid and data replacement arms, the across-grid install arm, the forced refresh arm and the data replacement with test rows arm, by their state checks |
  | the test rows' drop removed: in setData's two catches; in `setTestData` | a segfault in the data replacement with test rows arm; in the test data replacement arm |
  | the store's replacement hold removed | the data replacement arm's state check |
  | a column's grid installed in pieces; the derivation's hold removed | the forced refresh arm's grid check; the growing-grid arm's |
  | the rederive's settle removed | the linear rederive arm's state check |
  | the second round's engine mutants again (grow, prior tree and install marks, variance marks, map, table, monotone, vector fits, variance reset, compaction) | as in the table above, or by an earlier arm's state check or trap |

  Two survive. The variance refresh's per-tree mark removed: the mean forest's last mark still stands,
  so a throw there resets every variance tree, more than the mark would, and the chain is as
  consistent. The rebuild phase's in-place settle removed: no failure in any arm reaches that phase (a
  probe aborting in its catch never fires), which writes into storage the validate phase left sized, so
  the settle is a guard.

  The engine mutants of the second fix round were run against the alloc suite alone, built from the
  same sources; the state checks are the reads just after a failure, before any further sweep. The floor
  checks reuse one sampler for a run of edits, so once a mutant lets an edit in, the later
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
   gamma scale rounds to 0), at every door the floor reaches. Under a prior whose draw cannot reach the
   value it is refused, even in the own state of a sampler under another prior, since the floor is
   judged against the receiving chain's prior. Alternative: refuse such scales when the prior is built, a change to
   what a user may pass.
3. A warm-start donor whose DART weights are well formed but hold 0 on a column this forest may split on
   starts the recipient's DART fresh; setState refuses the same weights. Such a donor is package-written
   (a `forest(x1 + x2)` fit on the same three predictors, which the tinytest builds), and a warm start
   crosses model differences by design. Alternative: refuse it as `valueMismatch`.
4. After an allocation failure only the tree in flight is settled, rather than every constant-leaf
   tree's map (amendment 4's text): every other tree finished its move and draw, so its map and table
   agree. The settle makes the tree its state describes at once, its map rewritten from its partition,
   rather than marking the map stale for its next draw, so a state read just after a failure, or a
   sampler given it, has the failed chain's fits. A tree a failure stops mid-rebuild goes back to a
   single leaf, and the monotone tree in flight restarts its block at zero even when its values are in
   order, since checking the order allocates. Alternatives: mark the map stale, which leaves the state
   read before the next sweep describing another chain; mark every map stale; or keep a copy of a tree
   being rebuilt to put back, an allocation on every rebuild.
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
9. Amendment 16 left for a sweep: at any throw point the variance forest's combined variance equals
   the product of its factors within rounding, and the next sweep's closing recompute restores it
   exactly; a recompute in the catch would buy exactness only until that sweep. A variance tree a
   failure stops mid-draw or mid-install is different, its factors part written, and the settle
   recomputes the product there.
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
15. setState refuses a state whose DART weights hold 0 on a column this sampler's forest may split on,
    the same floor that refuses an edited 0, even when the state is another sampler's own, as a
    restricted-DART sampler's is when given to an unrestricted one. setState restores a state this
    sampler could have written, and none of its chains draws a weight of 0 on a column its forest may
    split on, the Dirichlet draw being positive there; the warm start, which crosses model differences,
    starts such weights fresh instead (call 3). Alternative: accept it, as main does, the next sweep's
    DART draw moving the weight off 0.
16. A move takes its own proposal back on a throw, from what it already holds (the birth's pair, the
    death's node, the subtree snapshot, the incumbent rule), and a birth's undo puts the structure
    back before it releases the pair, the free list's room reserved with the move's. Alternative: the
    settle collapses a proposal's empty leaves, which keeps a tree the chain never accepted.
17. A throw in a rebuild of every tree, in order, resets the tree in flight and every tree after it,
    the later forests' and the variance forest's included (`restUnreached_`), since they still index
    the rows or grid the rebuild replaced. Alternative: copy every tree first to put back, an
    allocation on every rebuild.
18. A chain the store has moved under while another chain threw, or before any chain was rebuilt,
    gets its trees reset (`Chain::resetEveryTree`) rather than re-routed: its rules still name the old
    grid's positions, and moving them allocates. Alternative: re-route without the move, which can
    leave a split past the new grid's last point.
19. In the unforced predictor transaction, a throw while the fits are rebuilt settles every tree in
    place, its structure and leaf values kept and its map and fits rewritten, every variance tree
    reset, rather than rolling the store back: by then the recovered parameters are partly written.
    Alternative: keep a copy of every leaf table and factor to roll back to.
20. The store's data replacement is all or nothing: its grids are derived over copies of the ones held
    and its codes written only after, at the cost of the grids' copy, never the codes'; a column's grid
    derived again or installed goes in whole the same way. Alternative: settle a half-replaced store, a
    column's count and points then apart.
21. A standardization rederive settles to the constants the leaf holds when the throw comes, rather
    than the ones before the call, which would need a copy made before anything moves. Alternative: that
    copy, an allocation on every rederive.
22. A data replacement with another row count and a rederive's saved draws are left unsettled (TODO,
    allocation-failure-residuals): settling the first means staging the response's and the leaf's
    per-row storage, every family's, before the store moves.
23. A throw while test rows are built or sized drops them, the store's test count back to 0, rather
    than keeping the rows held before: keeping them means building the new test store aside, a copy of
    every test column on every call. Alternative: that copy.

## Landing note

Landed by fast-forward at c202a6cf5 after four reviews, one over the build and one over each fix round, the third
round an audit of every path that rebuilds trees or moves the store. On the reviewed tip rebased onto main,
4f11d8d70: the nine per-push workflows by dispatch, eight green and lint failing on four mutation-battery anchors
(m02, m05, m26 and m27) that the moves' undo blocks had moved; re-anchored in c202a6cf5, each mutant killed again
at its gate, lint and doc-freshness green there. The x86 leg on 4f11d8d70 (Linux, x86-64): tests/cpp 548 checks
plain and under ASan/UBSan with the leak sanitizer, no report; the full tinytest suite, 42098 results, 0 failed;
the equivalence trio in statistical mode, 55 scenarios compared and indistinguishable, BCF and multinomial within
the cross-host bound. The last review's own gates on its reference build: the trio bitwise strict, 55, 19 and 11
identical; the snapshot files 14, 3, 7 and 3, 0 failed; tests/cpp under ASan/UBSan 548 ok. Two mutants survive,
both argued harmless (the variance refresh's per-tree mark, the rebuild phase's settle, a guard). Residuals in
TODO allocation-failure-residuals. The calls not put to the maintainer are dec-A266.

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
